#include <fstream>
#include <iostream>
#include <thread>
#include <queue>
#include <vector>
#include <errno.h>
#include <cstdlib>
#include <chrono>
#include <iomanip>

#include "TaskqHandler.h"
#include "SCIHarness.h"

#define MAX_SEND_CEV_SIZE 2000000

std::mutex print_duration_mutex;
std::mutex print_comm_mutex;

void print_duration(std::chrono::_V2::system_clock::time_point t1, string tag) {
    print_duration_mutex.lock();
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "::" << tag << " took "
              << ((double)std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()) / 1000
              << " seconds\n";
    print_duration_mutex.unlock();
}

void print_comm(uint64_t c_byte, string tag) {
    print_comm_mutex.lock();
    auto c_mbyte = ((double)c_byte) / (1024 * 1024);
    std::cout << ">>" << tag << " comm " << c_mbyte << " MB\n";
    print_comm_mutex.unlock();
}

void split_task_queue_on_destination(std::queue<Task>* taskq, std::vector<std::vector<Task>>& taskvs, std::vector<Task>& plaintext_taskv) {
	TaskComm& taskComm = TaskComm::getClientInstance();
	size_t tileNum = taskComm.getTileNum();
	size_t tileIndex = taskComm.getTileIndex();
    if (taskvs.size() != tileNum) taskvs.resize(tileNum);

	while (!taskq->empty()) {
		Task cur_task = taskq->front();
		bool isPlain = true;

        TaskPayload* tp = (TaskPayload*)cur_task.buf;
        for (int k=0; k<tp->operandNum; ++k)
            if (tp->isOperandEncrypted[k])
                isPlain = false;

		if (!isPlain) {
			taskvs[cur_task.dstTid].push_back(cur_task);
		} else {
			plaintext_taskv.push_back(cur_task);
		}

		taskq->pop();
	}

	// for (int i=0; i<taskvs.size(); i++) {
	// 	if (i != tileIndex) {
	// 		printf("taskvs[%d].size = %d\n", i, taskvs[i].size());
	// 	}
	// }

    // printf("plainv.size = %d\n", plaintext_taskv.size());

    // Move task digest from tileIndex to the encrypter
    for (Task task : taskvs[tileIndex]) {
        int cur_dstTid = -1;
        TaskPayload* tp = (TaskPayload*)task.buf;
        for (int k=0; k<tp->operandNum; ++k) {
            if (tp->isOperandEncrypted[k]) {
                cur_dstTid = tp->getEncTid(k);
                break;
            }
        }
        if (cur_dstTid != -1) {
            taskvs[cur_dstTid].push_back(task);
        }
    }	
    taskvs[tileIndex].clear();
}

void split_task_queue_on_destination(std::vector<Task>& taskq, std::vector<std::vector<Task>>& taskvs, std::vector<Task>& plaintext_taskv) {
	TaskComm& taskComm = TaskComm::getClientInstance();
	size_t tileNum = taskComm.getTileNum();
	size_t tileIndex = taskComm.getTileIndex();
    if (taskvs.size() != tileNum) taskvs.resize(tileNum);

	for (Task& cur_task : taskq) {
		bool isPlain = true;

        TaskPayload* tp = (TaskPayload*)cur_task.buf;
        for (int k=0; k<tp->operandNum; ++k)
            if (tp->isOperandEncrypted[k])
                isPlain = false;

		if (!isPlain) {
			taskvs[cur_task.dstTid].push_back(cur_task);
		} else {
			plaintext_taskv.push_back(cur_task);
		}
	}
    taskq.clear();

	// for (int i=0; i<taskvs.size(); i++) {
	// 	if (i != tileIndex) {
	// 		printf("taskvs[%d].size = %d\n", i, taskvs[i].size());
	// 	}
	// }

    // printf("plainv.size = %d\n", plaintext_taskv.size());

    // Move task digest from tileIndex to the encrypter
    for (Task task : taskvs[tileIndex]) {
        int cur_dstTid = -1;
        TaskPayload* tp = (TaskPayload*)task.buf;
        for (int k=0; k<tp->operandNum; ++k) {
            if (tp->isOperandEncrypted[k]) {
                cur_dstTid = tp->getEncTid(k);
                break;
            }
        }
        if (cur_dstTid != -1) {
            taskvs[cur_dstTid].push_back(task);
        }
    }	
    taskvs[tileIndex].clear();
}

void split_task_queue_on_destination(std::queue<Task>* taskq, std::vector<std::queue<Task>>& taskqs, std::queue<Task>& plaintext_taskq) {
	TaskComm& taskComm = TaskComm::getClientInstance();
	size_t tileNum = taskComm.getTileNum();
	size_t tileIndex = taskComm.getTileIndex();
    if (taskqs.size() != tileNum) taskqs.resize(tileNum);

	while (!taskq->empty()) {
		Task cur_task = taskq->front();
		bool isPlain = true;

        TaskPayload* tp = (TaskPayload*)cur_task.buf;
        for (int k=0; k<tp->operandNum; ++k)
            if (tp->isOperandEncrypted[k])
                isPlain = false;

		if (!isPlain) {
			taskqs[cur_task.dstTid].push(cur_task);
		} else {
			plaintext_taskq.push(cur_task);
		}

		taskq->pop();
	}

	for (int i=0; i<taskqs.size(); i++) {
		if (i != tileIndex) {
			printf("taskqs[%d].size = %d\n", i, taskqs[i].size());
		}
	}

    printf("plainq.size = %d\n", plaintext_taskq.size());

    // Move task digest from tileIndex to the encrypter
    while (!taskqs[tileIndex].empty()) {
        Task& task = taskqs[tileIndex].front();
        int cur_dstTid = -1;
        TaskPayload* tp = (TaskPayload*)task.buf;
        for (int k=0; k<tp->operandNum; ++k) {
            if (tp->isOperandEncrypted[k]) {
                cur_dstTid = tp->getEncTid(k);
                break;
            }
        }
        if (cur_dstTid != -1) {
            taskqs[cur_dstTid].push(task);
        }

        taskqs[tileIndex].pop();
    }
}

void task_queue_handler(std::queue<Task>* taskq) {
	TaskComm& taskComm = TaskComm::getClientInstance();
	size_t tileNum = taskComm.getTileNum();
	size_t tileIndex = taskComm.getTileIndex();

	CryptoUtil& cryptoUtil = CryptoUtil::getInstance();

    std::vector<std::vector<Task>>& taskvs = taskComm.getLocalTaskvs(); // Vector indexed by destination tile index
	std::vector<Task> plaintext_taskv; // Plain-text task vector

    split_task_queue_on_destination(taskq, taskvs, plaintext_taskv);

	// Send taskq digest
    for (int i=0; i<taskvs.size(); i++) {
        if (i != tileIndex && taskvs[i].size() != 0) {
            Task* empty_task = new Task();
            empty_task->type = taskvs[i][0].type;
            taskComm.sendTaskqDigest(*empty_task, taskvs[i].size(), i);
        }
        taskComm.setExpectedTaskNum(taskvs[i].size(), i);
    }

    // Pack cipher text
    // Note that the below four vectors are indexed by ce.tid
    std::vector<std::vector<CipherEntry*>> ceVecs;
    std::vector<std::vector<int32_t>> taskIdVecs;
    std::vector<std::vector<int32_t>> operandIdVecs;
    std::vector<std::vector<int32_t>> dstTidVecs;
    std::vector<std::vector<std::shared_ptr<PackedCipherEntry>>> pceVecs;
    ceVecs.resize(tileNum);
    taskIdVecs.resize(tileNum);
    operandIdVecs.resize(tileNum);
    dstTidVecs.resize(tileNum);
    pceVecs.resize(tileNum);
    int plainNumPerCE = 0;
    for (int i=0; i<taskvs.size(); i++) {
        // printf("---> send Task Num %d\n", taskvs[i].size());
        if (i != tileIndex) {
            // Inter edge task. The (cooperative) dstTid is the tile owning the dest vertex.
            int64_t cur_taskv_size = taskvs[i].size();
            for (int j= 0 ; j<cur_taskv_size; j++) {
                if (j % 100000 == 0) printf("Prepare pack share j = %d\n", j);
                auto& task = taskvs[i][j];
                TaskPayload* tp = (TaskPayload*)task.buf;
                plainNumPerCE = tp->plainNumPerCE;
                for (int k=0; k<tp->operandNum; ++k) {
                    if (tp->isOperandEncrypted[k]) {
                        CipherEntry* cur_ce_ptr = tp->getCipherEntryPtr(k);
                        size_t cur_ce_tid = cur_ce_ptr->tid;
                        ceVecs[cur_ce_tid].push_back(cur_ce_ptr);
                        taskIdVecs[cur_ce_tid].push_back(j);
                        operandIdVecs[cur_ce_tid].push_back(k);
                        dstTidVecs[cur_ce_tid].push_back(i);
                    } else {
                        taskComm.sendPlain(i, j, k);
                    }
                }
            }
        }
    }

    // printf("Begin plainNumPerCE %d\n", plainNumPerCE);

    for (int i=0; i<taskvs.size(); i++)
        if (i != tileIndex && ceVecs[i].size() != 0) {
            // printf("plainNumPerCE %d\n", plainNumPerCE);
            cryptoUtil.packCipherText(ceVecs[i], taskIdVecs[i], operandIdVecs[i], dstTidVecs[i], plainNumPerCE, pceVecs[i]);
        }

    // Split random share & Send packed ciphertext share decrypting
    for (int i=0; i<pceVecs.size(); i++) {
        if (i != tileIndex) {
            int64_t cur_pceVec_size = pceVecs[i].size();

            #pragma omp parallel for
            for (int j= 0 ; j<cur_pceVec_size; j++) {
                if (j % 10000 == 0) printf("Send share decrypting j = %d\n", j);
                cryptoUtil.splitRandomShareFromPackedCipherEntry(pceVecs[i][j]);
                // Write splitted share to the local task operands
                for (int m=0; m<pceVecs[i][j]->taskId.size(); ++m) {
                    auto& task = taskvs[i][pceVecs[i][j]->taskId[m]];
                    TaskPayload* tp = (TaskPayload*)task.buf;
                    int32_t cur_start = m * plainNumPerCE;
                    for (int kk=0; kk<plainNumPerCE; ++kk) {
                        tp->writeShareToOperand(pceVecs[i][j]->share[cur_start + kk], pceVecs[i][j]->operandId[m], kk);
                    }                
                }
                pceVecs[i][j]->share.clear();
                taskComm.sendDecryptPacked(pceVecs[i][j]);
            }
        }
    }

	for (int i=0; i<taskvs.size(); i++) {
		if (i != tileIndex && taskvs[i].size() > 0) {
			taskComm.sendEndDecrypt(i);
		}
	}

	// Multi-thread task queue computation using gc circuit
	std::vector<std::thread> threads;
	for (int i=0; i<taskvs.size(); i++) {
		if (i != tileIndex && taskvs[i].size() > 0) {
			threads.emplace_back([&taskvs, &taskComm, i] () {
                printf("Finish flag before computing %d\n", (int)taskComm.getLocalFinishqFlag(i));
				taskComm.computeTaskq(taskvs[i], taskComm.getTileIndex(), i, emp::ALICE);
			});
		}
	}

	std::thread plaintext_thread([&plaintext_taskv, &taskComm] () {
		taskComm.computePlaintextTaskq(plaintext_taskv);
	});

    std::vector<std::thread> recv_threads;
    issue_client_recv_threads(recv_threads);

    for (auto& thrd : recv_threads)
        thrd.join();

	for (auto& thrd : threads)
		thrd.join();

	plaintext_thread.join();
	
	// Encrypt share and merge share
	// Faster way for inter-edge update data pre-aggregation and inter-edge update data delivery

	// ## In algo kernel: update data transfer (send share, encrypt share, merge share)

	// Push task results
	for (int i=0; i<taskvs.size(); i++) {
		if (i != tileIndex) {
			for (auto task : taskvs[i]) {
				taskq->push(task);
			}
		}
	}

	for (int i=0; i<plaintext_taskv.size(); i++) {
		taskq->push(plaintext_taskv[i]);
	}

    taskComm.clean();

    plaintext_taskv.clear();
}

void single_destination_task_vector_handler(std::vector<Task>& taskv, size_t dstTid) {
	TaskComm& taskComm = TaskComm::getClientInstance();
    TaskqHandlerConfig& thc = taskComm.getTaskqHandlerConfig(dstTid);
	size_t tileNum = taskComm.getTileNum();
	size_t tileIndex = taskComm.getTileIndex();

	CryptoUtil& cryptoUtil = CryptoUtil::getInstance();

    std::vector<std::vector<Task>>& taskvs = taskComm.getLocalTaskvs(); // Vector indexed by destination tile index

    taskvs[dstTid].clear();
    std::swap(taskvs[dstTid], taskv);

	// Send taskq digest
    if (thc.sendTaskqDigest) {
        Task* empty_task = new Task();
        if (taskvs[dstTid].size() != 0) empty_task->type = taskvs[dstTid][0].type;
        taskComm.sendTaskqDigest(*empty_task, taskvs[dstTid].size(), dstTid);
    }

    taskComm.setExpectedTaskNum(taskvs[dstTid].size(), dstTid);

    if (taskvs[dstTid].size() == 0) {
        taskComm.clean(dstTid);
        return;
    }

    if (thc.sendOperand) {
        // Pack cipher text
        // Note that the below four vectors are indexed by ce.tid
        // Inter edge task. The (cooperative) dstTid is the tile owning the dest vertex.
        std::vector<CipherEntry*> ceVec;
        std::vector<int32_t> taskIdVec;
        std::vector<int32_t> operandIdVec;
        std::vector<int32_t> dstTidVec;
        std::vector<std::shared_ptr<PackedCipherEntry>> pceVec;
        int plainNumPerCE = 0;
        int64_t cur_taskv_size = taskvs[dstTid].size();
        for (int j= 0 ; j<cur_taskv_size; j++) {
            if (j % 100000 == 0) printf("Prepare pack share j = %d\n", j);
            auto& task = taskvs[dstTid][j];
            TaskPayload* tp = (TaskPayload*)task.buf;
            plainNumPerCE = tp->plainNumPerCE;
            for (int k=0; k<tp->operandNum; ++k) {
                if (tp->isOperandEncrypted[k]) {
                    CipherEntry* cur_ce_ptr = tp->getCipherEntryPtr(k);
                    size_t cur_ce_tid = cur_ce_ptr->tid;
                    ceVec.push_back(cur_ce_ptr);
                    taskIdVec.push_back(j);
                    operandIdVec.push_back(k);
                    dstTidVec.push_back(dstTid);
                } else {
                    taskComm.sendPlain(dstTid, j, k);
                }
            }
        }

        cryptoUtil.packCipherText(ceVec, taskIdVec, operandIdVec, dstTidVec, plainNumPerCE, pceVec);

        // Split random share & Send packed ciphertext share decrypting
        int64_t cur_pceVec_size = pceVec.size();

        #pragma omp parallel for
        for (int j= 0 ; j<cur_pceVec_size; j++) {
            if (j % 10000 == 0) printf("Send share decrypting j = %d\n", j);
            cryptoUtil.splitRandomShareFromPackedCipherEntry(pceVec[j]);
            // Write splitted share to the local task operands
            for (int m=0; m<pceVec[j]->taskId.size(); ++m) {
                auto& task = taskvs[dstTid][pceVec[j]->taskId[m]];
                TaskPayload* tp = (TaskPayload*)task.buf;
                int32_t cur_start = m * plainNumPerCE;
                for (int kk=0; kk<plainNumPerCE; ++kk) {
                    tp->writeShareToOperand(pceVec[j]->share[cur_start + kk], pceVec[j]->operandId[m], kk);
                }                
            }
            pceVec[j]->share.clear();
            taskComm.sendDecryptPacked(pceVec[j]);
        }

        taskComm.sendEndDecrypt(dstTid);
    }

    taskComm.doTaskOperandRotation(taskvs[dstTid], thc.rotation);

	std::thread compute_thread([&taskvs, &taskComm, &thc, dstTid] () {
        taskComm.computeTaskq(taskvs[dstTid], taskComm.getTileIndex(), dstTid, emp::ALICE, thc.rotation);
    });

    if (thc.mergeResult) {
        std::thread recv_thread([dstTid, &taskComm]() {
            int ret = 0;
            while (ret == 0 && !taskComm.getLocalFinishqFlag(dstTid)) {
                ret = taskComm.recvEntrance(dstTid);
            }
        });
        recv_thread.join();
        compute_thread.join();
    } else {
        compute_thread.join();
        for (Task& task : taskvs[dstTid]) {
            task.finished = true;
        }
    }

    std::swap(taskvs[dstTid], taskv);

    taskComm.clean(dstTid); // TODO
}

void plain_task_vector_handler(std::vector<Task>& plaintext_taskv) {
    TaskComm& taskComm = TaskComm::getClientInstance();
    taskComm.computePlaintextTaskq(plaintext_taskv);
}

void slot_extract(std::vector<Task>& dst_taskv, std::vector<Task>& src_taskv, int step) {
    if (step <= 0) return;
    for (int i=0 ; i<src_taskv.size(); i+=step) {
        dst_taskv.push_back(src_taskv[i]);
    }
}

uint64_t client_task_vector_direct_handler(std::vector<Task>& taskv, size_t dstTid) {
	TaskComm& taskComm = TaskComm::getClientInstance();
    TaskqHandlerConfig& thc = taskComm.getTaskqHandlerConfig(dstTid);
	size_t tileIndex = taskComm.getTileIndex();
    size_t tileNum = taskComm.getTileNum();
    uint32_t mpcBasePort = taskComm.getMPCBasePort();
    uint64_t comm = 0;

    std::vector<Task> new_taskv;

    if (!taskComm.getIsRotationBased()) { // round-halved
        if (taskv.size() !=0) {
            if (thc.rotation > 1) {
                slot_extract(new_taskv, taskv, thc.rotation);
            } else {
                new_taskv = taskv;
            }
            if (thc.rotation > 0) {
                taskComm.doTaskOperandRotation(new_taskv, 1);
                // taskComm.computeTaskq(new_taskv, tileIndex, dstTid, emp::ALICE, 1);
                comm = sci::computeTaskq(new_taskv, tileIndex, dstTid, mpcBasePort, tileNum, tileIndex, sci::ALICE, 1);
            } else {
                // taskComm.computeTaskq(new_taskv, tileIndex, dstTid, emp::ALICE, 0);
                comm = sci::computeTaskq(new_taskv, tileIndex, dstTid, mpcBasePort, tileNum, tileIndex, sci::ALICE, 0);
            }

            for (Task& task : taskv) {
                task.finished = true;
            }
        }
    } else {
        if (taskv.size() !=0) {
            new_taskv = taskv;
            if (thc.rotation > 0) {
                taskComm.doTaskOperandRotation(new_taskv, thc.rotation);
            }
            comm = sci::computeTaskq(new_taskv, tileIndex, dstTid, mpcBasePort, tileNum, tileIndex, sci::ALICE, thc.rotation);

            for (Task& task : taskv) {
                task.finished = true;
            }
        }        
    }

    return comm;
}

uint64_t server_task_vector_direct_handler(std::vector<Task>& taskv, size_t srcTid) {
	TaskComm& taskComm = TaskComm::getServerInstance();
    TaskqHandlerConfig& thc = taskComm.getTaskqHandlerConfig(srcTid);
	size_t tileIndex = taskComm.getTileIndex();
    size_t tileNum = taskComm.getTileNum();
    uint32_t mpcBasePort = taskComm.getMPCBasePort();
    uint64_t comm = 0;

    std::vector<Task> new_taskv;
    
    if (!taskComm.getIsRotationBased()) { // round-halved
        if (taskv.size() !=0) {
            if (thc.rotation > 1) {
                slot_extract(new_taskv, taskv, thc.rotation);
            } else {
                new_taskv = taskv;
            }
            if (thc.rotation > 0) {
                taskComm.doTaskOperandRotation(new_taskv, 1);
                // taskComm.computeTaskq(new_taskv, tileIndex, srcTid, emp::BOB, 1);
                comm = sci::computeTaskq(new_taskv, tileIndex, srcTid, mpcBasePort, tileNum, tileIndex, sci::BOB, 1);
            } else {
                // taskComm.computeTaskq(new_taskv, tileIndex, srcTid, emp::BOB, 0);
                comm = sci::computeTaskq(new_taskv, tileIndex, srcTid, mpcBasePort, tileNum, tileIndex, sci::BOB, 0);
            }

            for (Task& task : taskv) {
                task.finished = true;
            }
        }
    } else {
        if (taskv.size() !=0) {
            new_taskv = taskv;
            if (thc.rotation > 0) {
                taskComm.doTaskOperandRotation(new_taskv, thc.rotation);
            }
            comm = sci::computeTaskq(new_taskv, tileIndex, srcTid, mpcBasePort, tileNum, tileIndex, sci::BOB, thc.rotation);

            for (Task& task : taskv) {
                task.finished = true;
            }
        }        
    }

    return comm;
}

uint64_t get_rotation_upper_bound(std::vector<Task>& taskv) {
    // Assume that in taskv, tasks of the same destination vertex are grouped.
    uint64_t cur_group_vid = -1;
    uint64_t accumulated = 0;
    uint64_t max_accumulated = 0;
    for (Task& task : taskv) {
        if (task.vertexIndex == cur_group_vid) {
            accumulated += 1;
        } else {
            accumulated = 1;
            cur_group_vid = task.vertexIndex;
        }
        if (accumulated > max_accumulated) {
            max_accumulated = accumulated;
        }
    }

    uint64_t power = 1;
    while(power < max_accumulated)
        power <<= 1;
    return power;
}

uint64_t get_next_power_of_2(uint64_t num) {
    uint64_t power = 1;
    while(power < num)
        power <<= 1;
    return power;
}

void client_direct_recv(std::vector<Task>& taskv, size_t dstTid) {
	TaskComm& taskComm = TaskComm::getClientInstance();

    std::vector<std::vector<Task>>& taskvs = taskComm.getLocalTaskvs(); // Vector indexed by destination tile index

    taskvs[dstTid].clear();
    std::swap(taskvs[dstTid], taskv);

    taskComm.setExpectedTaskNum(taskvs[dstTid].size(), dstTid);

    if (taskvs[dstTid].size() == 0) {
        taskComm.clean(dstTid);
        return;
    }

    std::thread recv_thread([dstTid, &taskComm]() {
        int ret = 0;
        while (ret == 0 && !taskComm.getLocalFinishqFlag(dstTid)) {
            ret = taskComm.recvEntrance(dstTid);
        }
    });
    recv_thread.join();

    std::swap(taskvs[dstTid], taskv);

    taskComm.clean(dstTid);
}

void issue_client_recv_threads(std::vector<std::thread>& threads) {
	TaskComm& taskComm = TaskComm::getClientInstance();
	size_t tileNum = taskComm.getTileNum();
	size_t tileIndex = taskComm.getTileIndex();
    std::vector<std::vector<Task>>& taskvs = taskComm.getLocalTaskvs();

	for (int i = 0; i < tileNum; i++) {
		if (i != tileIndex && taskvs[i].size() > 0) {
			threads.emplace_back([i, &taskComm]() {
				int ret = 0;
                // int cnt = 0;
				while (ret == 0 && !taskComm.getLocalFinishqFlag(i)) {
					ret = taskComm.recvEntrance(i);
                    // cnt++;
				}
                // printf("client recv cnt %d\n", ret);
			});
		}
	}

	return;
}

void issue_server_recv_threads(std::vector<std::thread>& threads) {
	TaskComm& taskComm = TaskComm::getServerInstance();
	size_t tileNum = taskComm.getTileNum();
	size_t tileIndex = taskComm.getTileIndex();
    int party = emp::BOB;
    uint32_t mpcBasePort = taskComm.getMPCBasePort();

	for (int i = 0; i < tileNum; i++) {
		if (i != tileIndex) {
			threads.emplace_back([i, &taskComm, party, tileIndex, tileNum, mpcBasePort]() {
                uint32_t dstTid = i;
                // emp::NetIO * io = new emp::NetIO(party==emp::ALICE ? nullptr : "127.0.0.1", party==emp::ALICE? (mpcBasePort + dstTid*tileNum + tileIndex) : (mpcBasePort + tileIndex*tileNum + dstTid), false);
                // emp::setup_semi_honest(io, party);

				int ret = 0;
				while (ret == 0) {
					ret = taskComm.recvEntrance(i);
					if (taskComm.getRemoteFinishqFlag(i))
						taskComm.clean(i);
				}

                // emp::finalize_semi_honest();
                // delete io;
			});
		}
	}

	return;
}

void TaskComm::setUp(bool isClient) {
    isClient_ = isClient;

	basePort = 1812;
    if (!isCluster_) {
        ip = "127.0.0.1";
    } else {
        ip = "10.0.0.";
        ip += std::to_string(getTileIndex() + 1);
    }
    mpcBasePort = 12345;
    if (!isCluster_) {
        mpcIp = "127.0.0.1";
    } else {
        mpcIp = "10.0.0.";
        mpcIp += std::to_string(getTileIndex() + 1);
    }

    size_t tileNum = getTileNum();
    size_t tileIndex = getTileIndex();
    channels.resize(tileNum);

    std::vector<std::mutex> mtx_list(tileNum);
    channels_mtx.swap(mtx_list);
    std::vector<std::mutex> client_mtx_list(tileNum);
    client_computeTaskq_mtx.swap(client_mtx_list);    
    std::vector<std::condition_variable> cv_list(tileNum);
    client_computeTaskq_cv.swap(cv_list);
    std::vector<std::atomic<bool>> workdone_list(tileNum);
    client_workdone.swap(workdone_list);

    std::vector<Semaphore> client_smp_list(tileNum);
    client_computeTaskq_smp.swap(client_smp_list);

    std::vector<Semaphore> remote_update_ready_smp_list(tileNum);
    remote_update_ready_smp.swap(remote_update_ready_smp_list);

    std::vector<Semaphore> local_update_ready_smp_list(tileNum);
    local_update_ready_smp.swap(local_update_ready_smp_list);

    std::vector<std::atomic<uint64_t>> ready_tn_list(tileNum);
    readyTaskNum.swap(ready_tn_list); 

    std::vector<std::atomic<uint64_t>> finish_tn_list(tileNum);
    finishTaskNum.swap(finish_tn_list);
    
    taskvs.resize(tileNum);
    local_taskvs.resize(tileNum);
    expectedTaskNum.resize(tileNum);
    decrypt_taskq.resize(tileNum);
    merge_encrypted_taskq.resize(tileNum);
    packed_ce_q.resize(tileNum);

    remoteFinishqFlag.resize(tileNum);
    localFinishqFlag.resize(tileNum);

    thc_.resize(tileNum);
    
    for (int i=0; i<tileNum; ++i) {
        expectedTaskNum[i] = (uint64_t)-1;
        finishTaskNum[i] = 0;
        readyTaskNum[i] = 0;
        remoteFinishqFlag[i] = false;
        localFinishqFlag[i] = false;
        client_workdone[i] = false;
        // channels_mtx[i].unlock();
    }

    ioService = new osuCrypto::BtIOService(0);

    bool isHost = !isClient;

    // Establish channels
    for (uint32_t j = 0; j < tileNum; ++j) {
        if (j != tileIndex) {
            uint32_t port = 0;
            std::string name("endpoint:");
            if (isHost) {
                name += std::to_string(tileIndex) + "->" + std::to_string(j);
                port = basePort + (uint32_t)tileIndex;
                printf("Not host\n");
            }
            else
            {
                name += std::to_string(j) + "->" + std::to_string(tileIndex);
                port = basePort + (uint32_t)j;
            }
            printf("Endpoints name %s\n", name.c_str());
            if (!isCluster_ || isHost) {
                endpoints.emplace_back(*ioService, ip, port, isHost, name);
            } else {
                std::string remoteIp = "10.0.0.";
                remoteIp += std::to_string(j + 1);
                endpoints.emplace_back(*ioService, remoteIp, port, isHost, name);
            }
            channels[j] = &endpoints.back().addChannel("chl", "chl");
        }
    }

    if (isClient) {
        // Transfer public keys
        for (uint32_t j = 0; j < tileNum; ++j) {
            if (j != tileIndex) {
                PublicKeyEntry pke;
                memset(pke.pk, 0, MAX_PUBLIC_KEY_SIZE);
                cryptoUtil.serializePublicKeyToBytes(cryptoUtil.getPublicKey(), pke.pk);
                sendPublicKey(j, pke);
            }
        }
    } else {
        for (uint32_t j = 0; j < tileNum; ++j) {
            if (j != tileIndex) {
                recvEntrance(j);
            }
        }        
    }
}

void TaskComm::clean() {
    size_t tileNum = getTileNum();
    
    for (int i=0; i<tileNum; ++i) {
        expectedTaskNum[i] = (uint64_t)-1;
        finishTaskNum[i] = 0;
        readyTaskNum[i] = 0;
        remoteFinishqFlag[i] = false;
        localFinishqFlag[i] = false;
        // channels_mtx[i].unlock();
    } 

    for (int i=0; i<tileNum; ++i) {
        // for (auto& task : taskvs[i]) {
        //     task.delete_task_content_buf();
        // }
        taskvs[i].clear();
        // for (auto& task : local_taskvs[i]) {
        //     task.delete_task_content_buf();
        // }
        local_taskvs[i].clear();

        for (auto& task : decrypt_taskq[i]) {
            task.delete_task_content_buf();
        }
        decrypt_taskq[i].clear();

        for (auto& task : merge_encrypted_taskq[i]) {
            task.delete_task_content_buf();
        }
        merge_encrypted_taskq[i].clear();
    }
}

void TaskComm::clean(const uint32_t dstTid) {
    expectedTaskNum[dstTid] = (uint64_t)-1;
    finishTaskNum[dstTid] = 0;
    readyTaskNum[dstTid] = 0;
    remoteFinishqFlag[dstTid] = false;
    localFinishqFlag[dstTid] = false;

    for (auto& task : taskvs[dstTid]) {
        task.delete_task_content_buf();
    }
    taskvs[dstTid].clear();

    for (auto& task : decrypt_taskq[dstTid]) {
        task.delete_task_content_buf();
    }
    decrypt_taskq[dstTid].clear();

    for (auto& task : merge_encrypted_taskq[dstTid]) {
        task.delete_task_content_buf();
    }
    merge_encrypted_taskq[dstTid].clear();

    for (auto& task : local_taskvs[dstTid]) {
        task.delete_task_content_buf();
    }
    local_taskvs[dstTid].clear();
}

void TaskComm::closeChannels() {
    for (auto chl : channels)
        if (chl)
            chl->close();

    for (auto& ep : endpoints)
        ep.stop();

    ioService->stop();
}

void TaskComm::sendTaskDigest(Task task, const uint32_t dstTid) {
    std::string cmd_str = "TASK_DIGEST";

    std::string serial_str;
    boost::iostreams::back_insert_device<std::string> inserter(serial_str);
    boost::iostreams::stream<boost::iostreams::back_insert_device<std::string>> s(inserter);
    boost::archive::binary_oarchive oa(s);

    oa << cmd_str << task;

    s.flush();

    channels_mtx[dstTid].lock();
    channels[dstTid]->asyncSend(std::move(std::unique_ptr<osuCrypto::ByteStream>(new osuCrypto::ByteStream((uint8_t*)serial_str.data(), serial_str.size()))));
    channels_mtx[dstTid].unlock();
}

void TaskComm::sendTaskqDigest(Task task, const uint64_t taskNum, const uint32_t dstTid) {
    std::string cmd_str = "TASKQ_DIGEST";

    std::string serial_str;
    boost::iostreams::back_insert_device<std::string> inserter(serial_str);
    boost::iostreams::stream<boost::iostreams::back_insert_device<std::string>> s(inserter);
    boost::archive::binary_oarchive oa(s);

    oa << cmd_str << task << taskNum;

    s.flush();

    channels_mtx[dstTid].lock();
    channels[dstTid]->asyncSend(std::move(std::unique_ptr<osuCrypto::ByteStream>(new osuCrypto::ByteStream((uint8_t*)serial_str.data(), serial_str.size()))));
    channels_mtx[dstTid].unlock();    
}

void TaskComm::sendRotationUB(const uint64_t rotation_ub, const uint32_t dstTid) {
    std::string cmd_str = "ROTATION_UB";

    std::string serial_str;
    boost::iostreams::back_insert_device<std::string> inserter(serial_str);
    boost::iostreams::stream<boost::iostreams::back_insert_device<std::string>> s(inserter);
    boost::archive::binary_oarchive oa(s);

    oa << cmd_str << rotation_ub;

    s.flush();

    channels_mtx[dstTid].lock();
    channels[dstTid]->asyncSend(std::move(std::unique_ptr<osuCrypto::ByteStream>(new osuCrypto::ByteStream((uint8_t*)serial_str.data(), serial_str.size()))));
    channels_mtx[dstTid].unlock();    
}

void TaskComm::sendExpectedTaskNum(const uint64_t taskNum, const uint64_t dstTid) {
    std::string cmd_str = "EXPECTED_TASK_NUM";

    std::string serial_str;
    boost::iostreams::back_insert_device<std::string> inserter(serial_str);
    boost::iostreams::stream<boost::iostreams::back_insert_device<std::string>> s(inserter);
    boost::archive::binary_oarchive oa(s);

    oa << cmd_str << taskNum;

    s.flush();

    channels_mtx[dstTid].lock();
    channels[dstTid]->asyncSend(std::move(std::unique_ptr<osuCrypto::ByteStream>(new osuCrypto::ByteStream((uint8_t*)serial_str.data(), serial_str.size()))));
    channels_mtx[dstTid].unlock();
}

void TaskComm::taskDigestHandler(Task task, const uint32_t srcTid) {
    task.readyOperandNum = 0;
    switch (task.type) {
        case ADD_PAIR_DOUBLE_UINT:
            task.buf = (void*)(new AddPair<double, uint32_t>());
            break;
        case ADD_MIXED_PAIR_DOUBLE_UINT:
            task.buf = (void*)(new AddMixedPair<double, uint32_t>());
            break;
        case MIN_UINT_WITH_PARENT:
            task.buf = (void*)(new MinWithParent<uint32_t, uint64_t>());
            break;
        case ADD_DOUBLE:
            task.buf = (void*)(new Add<double>());
            break;
        case ADD_ULONG:
            task.buf = (void*)(new Add<uint64_t>());
            break;
        case ADD_UINT:
            task.buf = (void*)(new Add<uint32_t>());  
            break;   
        case DIV_DOUBLE:
            task.buf = (void*)(new Div<double>());
            break;
        case DEC_UINT: {
            task.buf = (void*)(new DecryptEntry<uint32_t>());
            break;
        }
        case ADD_UINT_WITH_REPLACE_PARENT: {
            task.buf = (void*)(new AddWithReplaceParent<uint32_t, uint64_t>());
            break;
        }
        case DEC_ULONG: {
            task.buf = (void*)(new DecryptEntry<uint64_t>());
            break;
        }
        case DEC_DOUBLE: {
            task.buf = (void*)(new DecryptEntry<double>());
            break;
        }
        case DEC_PAIR_UINT_ULONG: {
            task.buf = (void*)(new DecryptEntryPair<uint32_t, uint64_t>());
        }
        default:
            break;
    }

    taskvs[srcTid].push_back(task);
}

void TaskComm::taskqDigestHandler(Task emptyTask, const uint64_t taskNum, const uint32_t srcTid) {
    taskvs[srcTid].resize(taskNum);
    for (Task& task : taskvs[srcTid]) {
        task.readyOperandNum = 0;
        task.type = emptyTask.type;
        switch (task.type) {
            case ADD_PAIR_DOUBLE_UINT:
                task.buf = (void*)(new AddPair<double, uint32_t>());
                break;
            case ADD_MIXED_PAIR_DOUBLE_UINT:
                task.buf = (void*)(new AddMixedPair<double, uint32_t>());
                break;
            case MIN_UINT_WITH_PARENT:
                task.buf = (void*)(new MinWithParent<uint32_t, uint64_t>());
                break;
            case ADD_DOUBLE:
                task.buf = (void*)(new Add<double>());
                break;
            case ADD_ULONG:
                task.buf = (void*)(new Add<uint64_t>());
                break;
            case ADD_UINT:
                task.buf = (void*)(new Add<uint32_t>());  
                break;   
            case DIV_DOUBLE:
                task.buf = (void*)(new Div<double>());
                break;
            case DEC_UINT: {
                task.buf = (void*)(new DecryptEntry<uint32_t>());
                break;
            }
            case ADD_UINT_WITH_REPLACE_PARENT: {
                task.buf = (void*)(new AddWithReplaceParent<uint32_t, uint64_t>());
                break;
            }
            case DEC_ULONG: {
                task.buf = (void*)(new DecryptEntry<uint64_t>());
                break;
            }
            case DEC_DOUBLE: {
                task.buf = (void*)(new DecryptEntry<double>());
                break;
            }
            case DEC_PAIR_UINT_ULONG: {
                task.buf = (void*)(new DecryptEntryPair<uint32_t, uint64_t>());
            }
            default:
                break;
        }
    }
}

void TaskComm::sendDecrypt(CipherEntry& ce, const uint32_t dstTid, const uint32_t taskId, const uint32_t operandId0) {
    std::string cmd_str = "DECRYPT";

    std::string serial_str;
    boost::iostreams::back_insert_device<std::string> inserter(serial_str);
    boost::iostreams::stream<boost::iostreams::back_insert_device<std::string>> s(inserter);
    boost::archive::binary_oarchive oa(s);  

    oa << cmd_str << ce << dstTid << taskId << operandId0;

    s.flush();

    channels_mtx[ce.tid].lock();
    channels[ce.tid]->asyncSend(std::move(std::unique_ptr<osuCrypto::ByteStream>(new osuCrypto::ByteStream((uint8_t*)serial_str.data(), serial_str.size()))));   
    channels_mtx[ce.tid].unlock(); 
}

void TaskComm::sendDecryptPacked(std::shared_ptr<PackedCipherEntry> pce) {
    std::string cmd_str = "DECRYPT_PACKED";

    std::string serial_str;
    boost::iostreams::back_insert_device<std::string> inserter(serial_str);
    boost::iostreams::stream<boost::iostreams::back_insert_device<std::string>> s(inserter);
    boost::archive::binary_oarchive oa(s);  

    oa << cmd_str << *pce;

    s.flush();

    channels_mtx[pce->tid].lock();
    channels[pce->tid]->asyncSend(std::move(std::unique_ptr<osuCrypto::ByteStream>(new osuCrypto::ByteStream((uint8_t*)serial_str.data(), serial_str.size()))));   
    channels_mtx[pce->tid].unlock(); 
}

void TaskComm::sendEncrypted(CipherEntry& ce, const uint32_t dstTid, const uint32_t taskId, const uint32_t operandId0) {
    std::string cmd_str = "ENCRYPTED";

    std::string serial_str;
    boost::iostreams::back_insert_device<std::string> inserter(serial_str);
    boost::iostreams::stream<boost::iostreams::back_insert_device<std::string>> s(inserter);
    boost::archive::binary_oarchive oa(s);  

    oa << cmd_str << ce << taskId << operandId0;

    s.flush();

    channels_mtx[dstTid].lock();
    channels[dstTid]->asyncSend(std::move(std::unique_ptr<osuCrypto::ByteStream>(new osuCrypto::ByteStream((uint8_t*)serial_str.data(), serial_str.size()))));   
    channels_mtx[dstTid].unlock();   
}

void TaskComm::sendEndEncrypted(const uint32_t dstTid) {
    std::string cmd_str = "END_ENCRYPTED";

    std::string serial_str;
    boost::iostreams::back_insert_device<std::string> inserter(serial_str);
    boost::iostreams::stream<boost::iostreams::back_insert_device<std::string>> s(inserter);
    boost::archive::binary_oarchive oa(s);  

    oa << cmd_str;

    s.flush();

    channels_mtx[dstTid].lock();
    channels[dstTid]->asyncSend(std::move(std::unique_ptr<osuCrypto::ByteStream>(new osuCrypto::ByteStream((uint8_t*)serial_str.data(), serial_str.size()))));   
    channels_mtx[dstTid].unlock();   
}

void TaskComm::sendEndDecrypt(const uint32_t dstTid) {
    std::string cmd_str = "END_DECRYPT";

    std::string serial_str;
    boost::iostreams::back_insert_device<std::string> inserter(serial_str);
    boost::iostreams::stream<boost::iostreams::back_insert_device<std::string>> s(inserter);
    boost::archive::binary_oarchive oa(s);  

    oa << cmd_str;

    s.flush();

    channels_mtx[dstTid].lock();
    channels[dstTid]->asyncSend(std::move(std::unique_ptr<osuCrypto::ByteStream>(new osuCrypto::ByteStream((uint8_t*)serial_str.data(), serial_str.size()))));   
    channels_mtx[dstTid].unlock();   
}

void TaskComm::sendEndRotation(const uint32_t dstTid) {
    std::string cmd_str = "END_ROTATION";

    std::string serial_str;
    boost::iostreams::back_insert_device<std::string> inserter(serial_str);
    boost::iostreams::stream<boost::iostreams::back_insert_device<std::string>> s(inserter);
    boost::archive::binary_oarchive oa(s);  

    oa << cmd_str;

    s.flush();

    channels_mtx[dstTid].lock();
    channels[dstTid]->asyncSend(std::move(std::unique_ptr<osuCrypto::ByteStream>(new osuCrypto::ByteStream((uint8_t*)serial_str.data(), serial_str.size()))));   
    channels_mtx[dstTid].unlock();   
}

void TaskComm::sendPlain(const uint32_t dstTid, const uint32_t taskId, const uint32_t operandId0) {
    std::string cmd_str = "PLAIN";

    std::string serial_str;
    boost::iostreams::back_insert_device<std::string> inserter(serial_str);
    boost::iostreams::stream<boost::iostreams::back_insert_device<std::string>> s(inserter);
    boost::archive::binary_oarchive oa(s);  

    oa << cmd_str << taskId << operandId0;

    s.flush();

    channels_mtx[dstTid].lock();
    channels[dstTid]->asyncSend(std::move(std::unique_ptr<osuCrypto::ByteStream>(new osuCrypto::ByteStream((uint8_t*)serial_str.data(), serial_str.size()))));    
    channels_mtx[dstTid].unlock();
}

void TaskComm::decryptHandler(CipherEntry& ce, const uint32_t srcTid, const uint32_t dstTid, const uint32_t taskId, const uint32_t operandId0) {
    if (ce.tid != tileIndex_) {
        printf("Wrong send decrypt!\n");
        exit(-1);
    }
    // printf("ce.plainNum %d\n", ce.plainNum);
    if (ce.plainNum > 2 || ce.plainNum <= 0) {
        printf("Error plain num in CipherEntry!\n");
        exit(-1);
    }    
    
    ophelib::Ciphertext ct;
    cryptoUtil.deserializeCipherEntryToCipherText(ct, ce);
    ShareVec shareVec;
    cryptoUtil.decryptToInt64Vec(ct, shareVec.shares, ce.plainNum);
    // printf("decrypt share Vec %d, %lld, %lld\n", ce.plainNum, shareVec.shares[0], shareVec.shares[1]);

    if (dstTid == tileIndex_) {
        writeShareHandler(shareVec, srcTid, taskId, operandId0);
    } else {
        sendRouteShare(shareVec, srcTid, dstTid, taskId, operandId0);
    }

    return;
}

void TaskComm::decryptPackedHandler(const uint32_t srcTid) {
    #pragma omp parallel for
    for (int i=0; i<packed_ce_q[srcTid].size(); ++i) {
        std::shared_ptr<PackedCipherEntry>& pce = packed_ce_q[srcTid][i];
        if (pce->tid != tileIndex_) {
            printf("Wrong send packed decrypt!\n");
            exit(-1);
        }
        cryptoUtil.decryptPackedCipherEntry(pce);

        uint32_t ce_num = pce->taskId.size();
        uint32_t plainNumPerCE = pce->plainNum / ce_num;
        // printf("plainNumPerCE %d\n", plainNumPerCE);

        for (int m=0; m<pce->taskId.size(); ++m) {
            const uint32_t taskId = pce->taskId[m];
            const uint32_t dstTid = pce->dstTid[m];
            // printf("dstTid %d\n", dstTid);
            const uint32_t operandId = pce->operandId[m];
            int32_t cur_start = m * plainNumPerCE;
            ShareVec shareVec;
            for (int kk=0; kk<plainNumPerCE; ++kk) {
                shareVec.shares.push_back(pce->share[cur_start + kk]);
            }

            if (dstTid == tileIndex_) {
                writeShareHandler(shareVec, srcTid, taskId, operandId);
            } else {
                sendRouteShare(shareVec, srcTid, dstTid, taskId, operandId);
            }            
        }
    }

    packed_ce_q[srcTid].clear();
    
    return;
}

void TaskComm::mergeEncryptedHandler(CipherEntry& ce, const uint32_t srcTid, const uint32_t taskId, const uint32_t operandId0) {
    Task& task = local_taskvs[srcTid][taskId];
    uint8_t expectedOperandNum = 0;
    if (ce.tid != srcTid) {
        printf("Error encrypted share (encrypter) to be merged\n");
        exit(-1);
    }

    TaskPayload* tp = (TaskPayload*)task.buf;
    tp->mergeEncryptedShare(ce, operandId0);
    // tp->isOperandEncrypted[operandId0] = true;
    expectedOperandNum = 1;

    task.mtx.lock();
    task.readyOperandNum += 1;
    if (task.readyOperandNum.load() == expectedOperandNum) {
        task.finished = true;
        finishTaskNum[srcTid] += 1;
    }
    task.mtx.unlock();

    return;
}

void TaskComm::writeShareHandler(const ShareVec& shareVec, const uint32_t srcTid, const uint32_t taskId, const uint32_t operandId0) {
    Task& task = taskvs[srcTid][taskId];

    TaskPayload* tp = (TaskPayload*)task.buf;
    tp->isOperandEncrypted[operandId0] = true;
    for (int i=0; i<shareVec.shares.size(); ++i)
        tp->writeShareToOperand(shareVec.shares[i], operandId0, i);
    tp->setCipherEntryPlainNum(operandId0, shareVec.shares.size());
    uint8_t expectedOperandNum = tp->operandNum;

    task.mtx.lock();
    task.readyOperandNum += 1;
    if (task.readyOperandNum.load() == expectedOperandNum) readyTaskNum[srcTid] += 1;
    task.mtx.unlock();

    if (readyTaskNum[srcTid].load() == expectedTaskNum[srcTid]) {
		computeTaskq(getTaskv(srcTid), getTileIndex(), srcTid, emp::BOB);
        if (thc_[srcTid].mergeResult) sendEncryptedTaskqShare(srcTid);
        setRemoteFinishqFlag(srcTid);
    }

    return;
}

void TaskComm::writePlainHandler(const uint32_t srcTid, const uint32_t taskId, const uint32_t operandId0) {
    Task& task = taskvs[srcTid][taskId];
    // if (taskvs[srcTid].size() <= taskId) {
    //     printf("task.buf %lld\n", task.buf);
    //     printf("taskId %lld\n", taskId);
    //     printf("taskv length %lld\n", taskvs[srcTid].size());
    // }
    
    TaskPayload* tp = (TaskPayload*)task.buf;
    uint8_t expectedOperandNum = tp->operandNum;
    tp->isOperandEncrypted[operandId0] = false;

    task.mtx.lock();
    task.readyOperandNum += 1;
    if (task.readyOperandNum.load() == expectedOperandNum) readyTaskNum[srcTid] += 1;
    task.mtx.unlock();

    if (readyTaskNum[srcTid].load() == expectedTaskNum[srcTid]) {
		computeTaskq(getTaskv(srcTid), getTileIndex(), srcTid, emp::BOB);
        if (thc_[srcTid].mergeResult) sendEncryptedTaskqShare(srcTid);
        setRemoteFinishqFlag(srcTid);
    }

    return;
}

void TaskComm::sendRouteShare(const ShareVec& shareVec, const uint32_t srcTid, const uint32_t dstTid, const uint32_t taskId, const uint32_t operandId0) {
    std::string cmd_str = "ROUTE_SHARE";

    std::string serial_str;
    boost::iostreams::back_insert_device<std::string> inserter(serial_str);
    boost::iostreams::stream<boost::iostreams::back_insert_device<std::string>> s(inserter);
    boost::archive::binary_oarchive oa(s);  

    oa << cmd_str << shareVec << srcTid << taskId << operandId0;

    s.flush();

    channels_mtx[dstTid].lock();
    channels[dstTid]->asyncSend(std::move(std::unique_ptr<osuCrypto::ByteStream>(new osuCrypto::ByteStream((uint8_t*)serial_str.data(), serial_str.size())))); 
    channels_mtx[dstTid].unlock();      
}

int TaskComm::recvEntrance(const uint32_t srcTid) {
    osuCrypto::ByteStream recv;
    // channels_mtx[srcTid].lock();
    channels[srcTid]->recv(recv);
    // channels_mtx[srcTid].unlock();
    std::string serial_str((char*)recv.data(), recv.size());

    // wrap buffer inside a stream and deserialize serial_str into obj
    boost::iostreams::basic_array_source<char> device(serial_str.data(), serial_str.size());
    boost::iostreams::stream<boost::iostreams::basic_array_source<char> > s(device);
    boost::archive::binary_iarchive ia(s);

    std::string cmd_str;
    ia >> cmd_str;
    if (cmd_str == "TASK_DIGEST") {
        Task task;
        ia >> task;
        taskDigestHandler(task, srcTid);
    } else if (cmd_str == "ROTATION_UB") {
        TaskqHandlerConfig& thc = getTaskqHandlerConfig(srcTid);
        ia >> thc.rotation;
    } else if (cmd_str == "TASKQ_DIGEST") {
        Task task;
        uint64_t taskNum;
        ia >> task >> taskNum;
        if (taskNum == 0) {
            setRemoteFinishqFlag(srcTid);
        } else {
            taskqDigestHandler(task, taskNum, srcTid);
        }
    } else if (cmd_str == "DECRYPT") {
        CipherEntry ce;
        uint32_t dstTid;
        uint32_t taskId;
        uint32_t operandId0;
        ia >> ce >> dstTid >> taskId >> operandId0;
        Task task;
        task.buf = (void*)(new DecryptEntry<uint64_t>());
        task.srcTid = srcTid;
        task.dstTid = dstTid;
        task.vertexIndex = taskId; // tmp usage
        task.srcIndex = operandId0; // tmp usage
        DecryptEntry<uint64_t>* de = (DecryptEntry<uint64_t>*)task.buf;
        de->operand.enc = ce;
        decrypt_taskq[srcTid].push_back(task);
    } else if (cmd_str == "DECRYPT_PACKED") {
        std::shared_ptr<PackedCipherEntry> pce = std::make_shared<PackedCipherEntry>();
        ia >> *pce;
        packed_ce_q[srcTid].push_back(pce);
    } else if (cmd_str == "ROUTE_SHARE") {
        ShareVec shareVec;
        uint32_t srcTid;
        uint32_t taskId;
        uint32_t operandId0;
        ia >> shareVec >> srcTid >> taskId >> operandId0;
        writeShareHandler(shareVec, srcTid, taskId, operandId0);
    } else if (cmd_str == "EXPECTED_TASK_NUM") {
        uint64_t taskNum;
        ia >> taskNum;

        // Decrypt packed cipher entry
        decryptPackedHandler(srcTid);

        printf("received expected task num %d, received task num %d, ready tasknum %d\n", taskNum, taskvs[srcTid].size(), readyTaskNum[srcTid].load());
        
        expectedTaskNum[srcTid] = taskNum;
        
        if (readyTaskNum[srcTid].load() == expectedTaskNum[srcTid]) {
            computeTaskq(getTaskv(srcTid), getTileIndex(), srcTid, emp::BOB);
            if (thc_[srcTid].mergeResult) sendEncryptedTaskqShare(srcTid);
            setRemoteFinishqFlag(srcTid);
        }
    } else if (cmd_str == "END_DECRYPT") {
        decryptPackedHandler(srcTid);

        printf("received expected task num %d, ready tasknum %d\n", taskvs[srcTid].size(), readyTaskNum[srcTid].load());
        
        expectedTaskNum[srcTid] = taskvs[srcTid].size();
        
        if (readyTaskNum[srcTid].load() == expectedTaskNum[srcTid]) {
            computeTaskq(getTaskv(srcTid), getTileIndex(), srcTid, emp::BOB);
            if (thc_[srcTid].mergeResult) sendEncryptedTaskqShare(srcTid);
            setRemoteFinishqFlag(srcTid);
        }        
    } else if (cmd_str == "PLAIN") {
        uint32_t taskId;
        uint32_t operandId0;
        ia >> taskId >> operandId0;
        writePlainHandler(srcTid, taskId, operandId0);
    } else if (cmd_str == "ENCRYPTED") {
        CipherEntry ce;
        uint32_t taskId;
        uint32_t operandId0;
        ia >> ce >> taskId >> operandId0;

        Task task;
        task.buf = (void*)(new DecryptEntry<uint64_t>());
        task.srcTid = srcTid;
        task.vertexIndex = taskId; // tmp usage
        task.srcIndex = operandId0; // tmp usage
        DecryptEntry<uint64_t>* de = (DecryptEntry<uint64_t>*)task.buf;
        de->operand.enc = ce;
        merge_encrypted_taskq[srcTid].push_back(task);
        // mergeEncryptedHandler(ce, srcTid, taskId, operandId0);        
    } else if (cmd_str == "END_ENCRYPTED") {
        int64_t cur_taskv_size = this->merge_encrypted_taskq[srcTid].size();

        #pragma omp parallel for
        for (int j= 0 ; j<cur_taskv_size; j++) {
            if (j % 10000 == 0) printf("merging j = %d\n", j);
            auto& task = merge_encrypted_taskq[srcTid][j];
            DecryptEntry<uint64_t>* de = (DecryptEntry<uint64_t>*)task.buf;
            this->mergeEncryptedHandler(de->operand.enc, task.srcTid, task.vertexIndex, task.srcIndex);
        }

        for (auto& task : merge_encrypted_taskq[srcTid]) {
            task.delete_task_content_buf();
        }
        merge_encrypted_taskq[srcTid].clear();

        // printf("Here to set finish flag\n");
        if (finishTaskNum[srcTid] == expectedTaskNum[srcTid]) {
            setLocalFinishqFlag(srcTid);
        }

    } else if (cmd_str == "PUBLIC_KEY") {
        PublicKeyEntry pke;
        ia >> pke;
        ophelib::PublicKey cur_pk;
        cryptoUtil.deserializeBytesToPublicKey(cur_pk, pke.pk);
        cryptoUtil.setUpRemotePaillierCipher(srcTid, cur_pk);
    } else if (cmd_str == "END_ROTATION") {
        return -2;
    } else if (cmd_str == "FINISH") {
        return -1;
    }
        
    return 0;
}

void TaskComm::computeTaskq(vector<Task>& taskq, uint32_t curTid, uint32_t dstTid, int party, int rotation) {
    auto t1 = std::chrono::high_resolution_clock::now();

    printf("Begin compute taskq length %d %d\n", taskq.size(), party);
    printf("rotation %d\n", rotation);
    std::vector<std::thread> computeThreads;
    int compute_thread_num = 10;
    int64_t cur_taskv_size = taskq.size() - rotation;
    // print("cur_taskv_size %d\n")
    if (cur_taskv_size <= 0) {
        printf("End compute taskq length %d %d\n", taskq.size(), party);
        return;
    }

    int64_t step = (cur_taskv_size / compute_thread_num) + 1;
    
    for (int t=0; t<compute_thread_num; ++t) {
        computeThreads.emplace_back([&taskq, rotation, this, curTid, dstTid, t, step, cur_taskv_size, party, compute_thread_num] () {
            if (t*step >= cur_taskv_size) return;

            emp::NetIO * io = new emp::NetIO(party==emp::ALICE ? nullptr : "127.0.0.1", party==emp::ALICE? (mpcBasePort + (dstTid * compute_thread_num)*tileNum_ + compute_thread_num * tileIndex_ + t) : (mpcBasePort + (tileIndex_ * compute_thread_num)*tileNum_ + compute_thread_num * dstTid + t), true);
            emp::setup_semi_honest(io, party);
            // printf("inside Begin compute taskq length %d %d\n", taskq.size(), party);
            for (int j= t*step ; j<cur_taskv_size && j<(t+1)*step; j++) {
                Task& task = taskq[j];
                TaskPayload* tp = (TaskPayload*)task.buf;
                switch (task.type) {
                    case ADD_PAIR_DOUBLE_UINT: {
                        this->ssAddPairDoubleUint(task, party);
                    } break;
                    case ADD_MIXED_PAIR_DOUBLE_UINT: {
                        this->ssAddMixedPairDoubleUint(task, party);
                    } break;
                    case MIN_UINT_WITH_PARENT: {
                        this->ssMinUintWithParent(task, party);
                    } break;
                    case ADD_DOUBLE: {
                        this->ssAddDouble(task, party);
                    } break;
                    case ADD_ULONG: {
                        this->ssAddUlong(task, party);
                    } break;
                    case ADD_UINT: {
                        this->ssAddUint(task, party);
                    } break;
                    case ADD_UINT_WITH_REPLACE_PARENT: {
                        this->ssAddUintWithReplaceParent(task, party);
                    } break;
                    case DIV_DOUBLE: {
                        this->ssDivDouble(task, party);
                    } break;
                    case SWAP_CIPHER_ENTRY: {
                        this->ssSwapCipherEntry(task, party);
                    } break;
                    default: {
                    } break;
                }
                
                tp->isOperandEncrypted[0] = true;
                if (rotation != 0) {
                    tp->unifyOperand();
                }
            }

            emp::finalize_semi_honest();
            delete io;
        });
    }

    for (auto& thrd : computeThreads)
        thrd.join();

    printf("End compute taskq length %d %d\n", taskq.size(), party);

    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "::computetaskq took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() / 1000
              << " seconds\n";
    return;
}

void TaskComm::ssAddPairDoubleUint(Task& task, int party) {
	auto task_content = (AddPair<double, uint32_t>*)task.buf;

    struct {
		emp::Integer a; // Fixed point encoding for a
		emp::Integer b;
	} operands[3];

    // Convert share into emp-compatible type and merge shares
	for (int i=0; i<2; ++i) {
		if (!task_content->isOperandEncrypted[i]) {
            operands[i].a = cryptoUtil.doubleToEmpInteger(task_content->operands[i].plain_a, emp::ALICE);
            operands[i].b = cryptoUtil.uint32ToEmpInteger(task_content->operands[i].plain_b, emp::ALICE);
		} else {
            std::vector<int64_t> shareVec = {task_content->operands[i].enc.share[0], task_content->operands[i].enc.share[1]}; 
            emp::Integer packed_operand = cryptoUtil.shareVecToEmpInteger(shareVec, emp::ALICE) + cryptoUtil.shareVecToEmpInteger(shareVec, emp::BOB);
            operands[i].a = packed_operand >> 64;
            operands[i].b = packed_operand & emp::Integer(128, (uint64_t)0xffffffffffffffff, emp::PUBLIC);
 		}
	}

    // Compute
	operands[2].a = operands[0].a + operands[1].a;
	operands[2].b = operands[0].b + operands[1].b;

    // Mask operands[1]
    if (task_content->useMask) {
        emp::Bit maskBit(task_content->operandMask && (!task.isDummy), emp::ALICE);
        operands[0].a = operands[0].a.select(maskBit, operands[2].a);
        operands[0].b = operands[0].b.select(maskBit, operands[2].b);
    } else {
        operands[0].a = operands[2].a;
        operands[0].b = operands[2].b;
    }
	
	// Generate new shares
    int64_t new_share_a0 = cryptoUtil.getRandomShare();
    emp::Integer share_a0 = cryptoUtil.shareToEmpInteger(new_share_a0, emp::ALICE);
    emp::Integer share_a1 = operands[0].a - share_a0;
    int64_t new_share_a1 = cryptoUtil.empIntegerToShare(share_a1, emp::BOB);

    int64_t new_share_b0 = cryptoUtil.getRandomShare();
    emp::Integer share_b0 = cryptoUtil.shareToEmpInteger(new_share_b0, emp::ALICE);
    emp::Integer share_b1 = operands[0].b - share_b0;
    int64_t new_share_b1 = cryptoUtil.empIntegerToShare(share_b1, emp::BOB);
	
    // Store new shares
    if (party == emp::ALICE) {
        task_content->operands[0].enc.share[0] = new_share_a0;
        task_content->operands[0].enc.share[1] = new_share_b0;
    }

    if (party == emp::BOB) {
        task_content->operands[0].enc.share[0] = new_share_a1;
        task_content->operands[0].enc.share[1] = new_share_b1;
    }

    return;
}

void TaskComm::ssAddMixedPairDoubleUint(Task& task, int party) {
	auto task_content = (AddMixedPair<double, uint32_t>*)task.buf;

    struct {
		emp::Float a; // Fixed point encoding for a
		uint32_t b;
	} operands[3];

    // Convert share into emp-compatible type and merge shares
	for (int i=0; i<2; ++i) {
		if (!task_content->isOperandEncrypted[i]) {
            operands[i].a = cryptoUtil.shareToEmpFloat(task_content->operands[i].plain_a, emp::ALICE);
		} else {
            operands[i].a = cryptoUtil.shareToEmpFloat(task_content->operands[i].enc_a.share[0], emp::ALICE) + cryptoUtil.shareToEmpFloat(task_content->operands[i].enc_a.share[0], emp::BOB);
 		}
        operands[i].b = task_content->operands[i].plain_b;
	}

    // Compute
	operands[2].a = operands[0].a + operands[1].a;
	operands[2].b = operands[0].b + operands[1].b;

    // Mask operands[1]
    if (task_content->useMask) {
        emp::Bit maskBit(task_content->operandMask && (!task.isDummy), emp::ALICE);
        operands[0].a = operands[0].a.If(maskBit, operands[2].a);
        if (task_content->operandMask && (!task.isDummy)) {
            operands[0].b = operands[2].b;
        }
    } else {
        operands[0].a = operands[2].a;
        operands[0].b = operands[2].b;
    }
	
	// Generate new shares
    int64_t new_share_a0 = cryptoUtil.getRandomShare();
    emp::Float share_a0 = cryptoUtil.shareToEmpFloat(new_share_a0, emp::ALICE);
    emp::Float share_a1 = operands[0].a - share_a0;
    int64_t new_share_a1 = cryptoUtil.empFloatToShare(share_a1, emp::BOB);
	
    // Store new shares
    if (party == emp::ALICE) {
        task_content->operands[0].enc_a.share[0] = new_share_a0;
    }

    if (party == emp::BOB) {
        task_content->operands[0].enc_a.share[0] = new_share_a1;
    }

    task_content->operands[0].plain_b = operands[0].b;

    return;
}

void TaskComm::ssMinUintWithParent(Task& task, int party) {
	auto task_content = (MinWithParent<uint32_t, uint64_t>*)task.buf;

	struct {
		emp::Integer val;
		emp::Integer parent;
	} operands[3];

	for (int i=0; i<2; ++i) {
		if (!task_content->isOperandEncrypted[i]) {
            operands[i].val = emp::Integer(64, task_content->operands[i].plain_val, emp::ALICE);
            operands[i].parent = emp::Integer(64, task_content->operands[i].plain_parent, emp::ALICE);
            // printf("Illegal unecrypted\n");
		} else {
            operands[i].val = emp::Integer(64, task_content->operands[i].enc.share[0], emp::ALICE) + emp::Integer(64, task_content->operands[i].enc.share[0], emp::BOB);
            operands[i].parent = emp::Integer(64, task_content->operands[i].enc.share[1], emp::ALICE) + emp::Integer(64, task_content->operands[i].enc.share[1], emp::BOB);
 		}
	}

    // int64_t a1 = cryptoUtil.empIntegerToShare(operands[0].val, emp::PUBLIC);
    // int64_t a2 = cryptoUtil.empIntegerToShare(operands[1].val, emp::PUBLIC);
    // int64_t a3 = cryptoUtil.empIntegerToShare(operands[0].parent, emp::PUBLIC);
    // int64_t a4 = cryptoUtil.empIntegerToShare(operands[1].parent, emp::PUBLIC);
    

    // if (party == emp::ALICE) {
    //     printf("Compute SS MWP\n");
    //     printf("vertex %lld\n", task.vertexIndex);
    //     printf("SS val[0] %lld\n", a1);
    //     printf("SS val[1] %lld\n", a2);
    //     printf("SS parent[0] %lld\n", a3);
    //     printf("SS parent[1] %lld\n", a4);
    // }

    // Compute
    emp::Bit selectBit = operands[0].val.geq(operands[1].val);
	operands[2].val = operands[0].val.select(selectBit, operands[1].val);
	operands[2].parent = operands[0].parent.select(selectBit, operands[1].parent);	

    // Mask operands[1]
    if (task_content->useMask) {
        // std::cout<<task_content->operandMask<<"++++"<<task.isDummy<<std::endl;
        emp::Bit maskBit(task_content->operandMask && (!task.isDummy), emp::ALICE);
        //emp::Bit maskBit(true, emp::ALICE);
        operands[0].val = operands[0].val.select(maskBit, operands[2].val);
        operands[0].parent = operands[0].parent.select(maskBit, operands[2].parent);
    } else {
        operands[0].val = operands[2].val;
        operands[0].parent = operands[2].parent;
    }

    // a1 = cryptoUtil.empIntegerToShare(operands[0].val, emp::PUBLIC);
    // a3 = cryptoUtil.empIntegerToShare(operands[0].parent, emp::PUBLIC);

    // if (party == emp::ALICE) {
    //     printf("Chosen\n");
    //     printf("SS val[0] %lld\n", a1);
    //     printf("SS parent[0] %lld\n", a3);
    // }

	// Generate new shares
    int64_t new_share_val0 = cryptoUtil.getRandomShare();
    int64_t new_share_parent0 = cryptoUtil.getRandomShare();
    emp::Integer share_val0 = cryptoUtil.shareToEmpInteger(new_share_val0, emp::ALICE);
    emp::Integer share_parent0 = cryptoUtil.shareToEmpInteger(new_share_parent0, emp::ALICE);

    emp::Integer share_val1 = operands[0].val - share_val0;
    emp::Integer share_parent1 = operands[0].parent - share_parent0;  

    int64_t new_share_val1 = cryptoUtil.empIntegerToShare(share_val1, emp::BOB);    
    int64_t new_share_parent1 = cryptoUtil.empIntegerToShare(share_parent1, emp::BOB);
	
    // Store new shares
    if (party == emp::ALICE) {
        task_content->operands[0].enc.share[0] = new_share_val0;
        task_content->operands[0].enc.share[1] = new_share_parent0;
    }
    if (party == emp::BOB) {
        task_content->operands[0].enc.share[0] = new_share_val1;
        task_content->operands[0].enc.share[1] = new_share_parent1;
    }
}

void TaskComm::ssAddUintWithReplaceParent(Task& task, int party) {
	auto task_content = (AddWithReplaceParent<uint32_t, uint64_t>*)task.buf;

	struct {
		emp::Integer val;
		emp::Integer parent;
	} operands[3];

	for (int i=0; i<2; ++i) {
		if (!task_content->isOperandEncrypted[i]) {
            operands[i].val = emp::Integer(64, task_content->operands[i].plain_val, emp::ALICE);
            operands[i].parent = emp::Integer(64, task_content->operands[i].plain_parent, emp::ALICE);
		} else {
            operands[i].val = emp::Integer(64, task_content->operands[i].enc.share[0], emp::ALICE) + emp::Integer(64, task_content->operands[i].enc.share[0], emp::BOB);
            operands[i].parent = emp::Integer(64, task_content->operands[i].enc.share[1], emp::ALICE) + emp::Integer(64, task_content->operands[i].enc.share[1], emp::BOB);
 		}
	}

    // printf("SS val[0] %lld\n", cryptoUtil.empIntegerToShare(operands[0].val, emp::PUBLIC));

    // Compute
	operands[2].val = operands[0].val + operands[1].val;
	operands[2].parent = operands[1].parent;

    // Mask operands[1]
    if (task_content->useMask) {
        emp::Bit maskBit(task_content->operandMask && (!task.isDummy), emp::ALICE);
        operands[0].val = operands[0].val.select(maskBit, operands[2].val);
        operands[0].parent = operands[0].parent.select(maskBit, operands[2].parent);
    } else {
        operands[0].val = operands[2].val;
        operands[0].parent = operands[2].parent;
    }

    // printf("Compute SS AWRP\n");
    // printf("SS val[0] %lld\n", cryptoUtil.empIntegerToShare(operands[0].val, emp::PUBLIC));
    // printf("SS parent[0] %lld\n", cryptoUtil.empIntegerToShare(operands[0].parent, emp::PUBLIC));

	// Generate new shares
    int64_t new_share_val0 = cryptoUtil.getRandomShare();
    int64_t new_share_parent0 = cryptoUtil.getRandomShare();
    emp::Integer share_val0 = cryptoUtil.shareToEmpInteger(new_share_val0, emp::ALICE);
    emp::Integer share_parent0 = cryptoUtil.shareToEmpInteger(new_share_parent0, emp::ALICE);

    emp::Integer share_val1 = operands[0].val - share_val0;
    emp::Integer share_parent1 = operands[0].parent - share_parent0;  

    int64_t new_share_val1 = cryptoUtil.empIntegerToShare(share_val1, emp::BOB);    
    int64_t new_share_parent1 = cryptoUtil.empIntegerToShare(share_parent1, emp::BOB);
	
    // Store new shares
    if (party == emp::ALICE) {
        task_content->operands[0].enc.share[0] = new_share_val0;
        task_content->operands[0].enc.share[1] = new_share_parent0;
    }
    if (party == emp::BOB) {
        task_content->operands[0].enc.share[0] = new_share_val1;
        task_content->operands[0].enc.share[1] = new_share_parent1;
    }
}



void TaskComm::ssAddDouble(Task& task, int party) {
	auto task_content = (Add<double>*)task.buf;

    emp::Float operands[3];

    // Convert share into emp-compatible type and merge shares
	for (int i=0; i<2; ++i) {
		if (!task_content->isOperandEncrypted[i]) {
			operands[i] = emp::Float(task_content->operands[i].plain, emp::ALICE);
		} else {
            operands[i] = cryptoUtil.shareToEmpFloat(task_content->operands[i].enc.share[0], emp::ALICE) + cryptoUtil.shareToEmpFloat(task_content->operands[i].enc.share[0], emp::BOB);
		}
	}

    // Compute
	operands[2] = operands[0] + operands[1];

    // Mask operands[1]
    if (task_content->useMask) {
        emp::Bit maskBit(task_content->operandMask && (!task.isDummy), emp::ALICE);
        operands[0] = operands[0].If(maskBit, operands[2]);
    } else {
        operands[0] = operands[2];
    }
	
	// Generate new shares
    int64_t new_share_0 = cryptoUtil.getRandomShare();
    emp::Float share_0 = cryptoUtil.shareToEmpFloat(new_share_0, emp::ALICE);
    emp::Float share_1 = operands[0] - share_0;
    int64_t new_share_1 = cryptoUtil.empFloatToShare(share_1, emp::BOB);
	
    // Store new shares
    if (party == emp::ALICE)
        task_content->operands[0].enc.share[0] = new_share_0;
    if (party == emp::BOB)
        task_content->operands[0].enc.share[0] = new_share_1;

    return;
}

void TaskComm::ssAddUlong(Task& task, int party) {
	auto task_content = (Add<uint64_t>*)task.buf;

    emp::Integer operands[3];

    // Convert share into emp-compatible type and merge shares
	for (int i=0; i<2; ++i) {
		if (!task_content->isOperandEncrypted[i]) {
			operands[i] = emp::Integer(64, task_content->operands[i].plain, emp::ALICE);
		} else {
            operands[i] = cryptoUtil.shareToEmpInteger(task_content->operands[i].enc.share[0], emp::ALICE) + cryptoUtil.shareToEmpInteger(task_content->operands[i].enc.share[0], emp::BOB);
		}
	}

    // Compute
	operands[2] = operands[0] + operands[1];

    // Mask operands[1]
    if (task_content->useMask) {
        emp::Bit maskBit(task_content->operandMask && (!task.isDummy), emp::ALICE);
        operands[0] = operands[0].select(maskBit, operands[2]);
    } else {
        operands[0] = operands[2];
    }
	
	// Generate new shares
    int64_t new_share_0 = cryptoUtil.getRandomShare();
    emp::Integer share_0 = cryptoUtil.shareToEmpInteger(new_share_0, emp::ALICE);
    emp::Integer share_1 = operands[0] - share_0;
    int64_t new_share_1 = cryptoUtil.empIntegerToShare(share_1, emp::BOB);
	
    // Store new shares
    if (party == emp::ALICE)
        task_content->operands[0].enc.share[0] = new_share_0;
    if (party == emp::BOB)
        task_content->operands[0].enc.share[0] = new_share_1;

    return;
}

void TaskComm::ssAddUint(Task& task, int party) {
	auto task_content = (Add<uint32_t>*)task.buf;

    emp::Integer operands[3];

    // Convert share into emp-compatible type and merge shares
	for (int i=0; i<2; ++i) {
		if (!task_content->isOperandEncrypted[i]) {
			operands[i] = emp::Integer(64, task_content->operands[i].plain, emp::ALICE);
		} else {
            operands[i] = cryptoUtil.shareToEmpInteger(task_content->operands[i].enc.share[0], emp::ALICE) + cryptoUtil.shareToEmpInteger(task_content->operands[i].enc.share[0], emp::BOB);
		}
	}

    // Compute
	operands[2] = operands[0] + operands[1];

    // Mask operands[1]
    if (task_content->useMask) {
        emp::Bit maskBit(task_content->operandMask && (!task.isDummy), emp::ALICE);
        operands[0] = operands[0].select(maskBit, operands[2]);
    } else {
        operands[0] = operands[2];
    }
	
	// Generate new shares
    int64_t new_share_0 = cryptoUtil.getRandomShare();
    emp::Integer share_0 = cryptoUtil.shareToEmpInteger(new_share_0, emp::ALICE);
    emp::Integer share_1 = operands[0] - share_0;
    int64_t new_share_1 = cryptoUtil.empIntegerToShare(share_1, emp::BOB);
	
    // Store new shares
    if (party == emp::ALICE)
        task_content->operands[0].enc.share[0] = new_share_0;
    if (party == emp::BOB)
        task_content->operands[0].enc.share[0] = new_share_1;

    return;
}

void TaskComm::ssDivDouble(Task& task, int party) {
	auto task_content = (struct Add<double>*)task.buf;

    emp::Float operands[3];

    // Convert share into emp-compatible type and merge shares
	for (int i=0; i<2; ++i) {
		if (!task_content->isOperandEncrypted[i]) {
			operands[i] = emp::Float(task_content->operands[i].plain, emp::ALICE);
		} else {
            operands[i] = cryptoUtil.shareToEmpFloat(task_content->operands[i].enc.share[0], emp::ALICE) + cryptoUtil.shareToEmpFloat(task_content->operands[i].enc.share[0], emp::BOB);
		}
	}

    // Compute
	operands[2] = operands[0] / operands[1];

    // Mask operands[1]
    if (task_content->useMask) {
        emp::Bit maskBit(task_content->operandMask && (!task.isDummy), emp::ALICE);
        operands[0] = operands[0].If(maskBit, operands[2]);
    } else {
        operands[0] = operands[2];
    }
	
	// Generate new shares
    int64_t new_share_0 = cryptoUtil.getRandomShare();
    emp::Float share_0 = cryptoUtil.shareToEmpFloat(new_share_0, emp::ALICE);
    emp::Float share_1 = operands[0] - share_0;
    int64_t new_share_1 = cryptoUtil.empFloatToShare(share_1, emp::BOB);
	
    // Store new shares
    if (party == emp::ALICE)
        task_content->operands[0].enc.share[0] = new_share_0;
    if (party == emp::BOB)
        task_content->operands[0].enc.share[0] = new_share_1;

    return;
}

void TaskComm::ssSwapCipherEntry(Task& task, int party) {
	auto task_content = (SwapCipherEntry*)task.buf;

	struct {
		emp::Integer val0;
		emp::Integer val1;
	} operands[3];

	for (int i=0; i<2; ++i) {
        operands[i].val0 = emp::Integer(64, task_content->operands[i].enc.share[0], emp::ALICE) + emp::Integer(64, task_content->operands[i].enc.share[0], emp::BOB);
        operands[i].val1 = emp::Integer(64, task_content->operands[i].enc.share[1], emp::ALICE) + emp::Integer(64, task_content->operands[i].enc.share[1], emp::BOB);
	}

    // Compute
    operands[2].val0 = operands[0].val0;
    operands[2].val1 = operands[0].val1;
    emp::Bit swapBit(task_content->operandMask, emp::ALICE);
	operands[0].val0 = operands[0].val0.select(swapBit, operands[1].val0);
    operands[0].val1 = operands[0].val1.select(swapBit, operands[1].val1);
	operands[1].val0 = operands[1].val0.select(swapBit, operands[2].val0);
    operands[1].val1 = operands[1].val1.select(swapBit, operands[2].val1);

	// Generate new shares
    int64_t new_share_val00 = cryptoUtil.getRandomShare();
    int64_t new_share_val10 = cryptoUtil.getRandomShare();
    emp::Integer share_val00 = cryptoUtil.shareToEmpInteger(new_share_val00, emp::ALICE);
    emp::Integer share_val10 = cryptoUtil.shareToEmpInteger(new_share_val10, emp::ALICE);

    emp::Integer share_val01 = operands[0].val0 - share_val00;
    emp::Integer share_val11 = operands[0].val1 - share_val10;  

    int64_t new_share_val01 = cryptoUtil.empIntegerToShare(share_val01, emp::BOB);    
    int64_t new_share_val11 = cryptoUtil.empIntegerToShare(share_val11, emp::BOB);
	
    // Store new shares
    if (party == emp::ALICE) {
        task_content->operands[0].enc.share[0] = new_share_val00;
        task_content->operands[0].enc.share[1] = new_share_val10;
    }
    if (party == emp::BOB) {
        task_content->operands[0].enc.share[0] = new_share_val01;
        task_content->operands[0].enc.share[1] = new_share_val11;
    }
}

void TaskComm::computePlaintextTaskq(vector<Task>& taskq) {
    for (Task& task : taskq) {
        switch (task.type) {
            case ADD_PAIR_DOUBLE_UINT: {
                AddPair<double, uint32_t>* addp = (AddPair<double, uint32_t>*)task.buf;
                addp->operands[0].plain_a += addp->operands[1].plain_a;
                addp->operands[0].plain_b += addp->operands[1].plain_b;
                break;
            }
            case ADD_MIXED_PAIR_DOUBLE_UINT: {
                AddMixedPair<double, uint32_t>* addmp = (AddMixedPair<double, uint32_t>*)task.buf;
                addmp->operands[0].plain_a += addmp->operands[1].plain_a;
                addmp->operands[0].plain_b += addmp->operands[1].plain_b;
                break;
            }
            case MIN_UINT_WITH_PARENT: {
                MinWithParent<uint32_t, uint64_t>* mwp = (MinWithParent<uint32_t, uint64_t>*)task.buf;
                if (mwp->operands[0].plain_val > mwp->operands[1].plain_val) {
                    mwp->operands[0].plain_val = mwp->operands[1].plain_val;
                    mwp->operands[0].plain_parent = mwp->operands[1].plain_parent;
                }
                break;
            }
            case ADD_DOUBLE: {
                Add<double>* add = (Add<double>*)task.buf;
                add->operands[0].plain += add->operands[1].plain;
                break;
            }
            case ADD_ULONG:{
                Add<uint64_t>* add = (Add<uint64_t>*)task.buf;
                add->operands[0].plain += add->operands[1].plain;
                break;
            }
            case ADD_UINT: {
                Add<uint32_t>* add = (Add<uint32_t>*)task.buf;
                add->operands[0].plain += add->operands[1].plain;
                break;
            } 
            case ADD_UINT_WITH_REPLACE_PARENT: {
                AddWithReplaceParent<uint32_t, uint64_t>* awrp = (AddWithReplaceParent<uint32_t, uint64_t>*)task.buf;
                awrp->operands[0].plain_val += awrp->operands[1].plain_val;
                awrp->operands[0].plain_parent = awrp->operands[1].plain_parent;
                break;
            } 
            case DIV_DOUBLE: {
                struct Div<double>* div = (struct Div<double>*)task.buf;
                div->operands[0].plain /= div->operands[1].plain;
                break;
            }
            default:
                break;
        }
		task.finished = true;
    }

    return;
}

void TaskComm::sendEncryptedTaskqShare(const uint32_t dstTid) {
    // printf("Begin send encrypted taskq share\n");
    #pragma omp parallel for
    for (int j=0; j<taskvs[dstTid].size(); ++j) {
        Task task = taskvs[dstTid][j];
        TaskPayload* tp = (TaskPayload*)task.buf;
        CipherEntry cur_ce = tp->encryptShare(0, tileIndex_);
        sendEncrypted(cur_ce, dstTid, j, 0);
    }

    sendEndEncrypted(dstTid);    

    return;
}

void TaskComm::sendEncryptedTaskqShare(const std::vector<Task>& taskv, const uint32_t dstTid) {
    #pragma omp parallel for
    for (int j=0; j<taskv.size(); ++j) {
        Task task = taskv[j];
        TaskPayload* tp = (TaskPayload*)task.buf;
        CipherEntry cur_ce = tp->encryptShare(0, tileIndex_);
        sendEncrypted(cur_ce, dstTid, j, 0);
    }

    sendEndEncrypted(dstTid);    

    return;
}

void TaskComm::setRemoteFinishqFlag(const uint32_t dstTid) {
    remoteFinishqFlag[dstTid] = true;
}

bool TaskComm::getRemoteFinishqFlag(const uint32_t dstTid) {
    return remoteFinishqFlag[dstTid];
}

void TaskComm::setLocalFinishqFlag(const uint32_t dstTid) {
    localFinishqFlag[dstTid] = true;
}

bool TaskComm::getLocalFinishqFlag(const uint32_t dstTid) {
    return localFinishqFlag[dstTid];
}

std::vector<std::thread>& TaskComm::getComputeRemoteTaskqThreads() {
    return computeRemoteTaskqThreads;
}

void TaskComm::sendPublicKey(const uint32_t dstTid, const PublicKeyEntry& pke) {
    std::string cmd_str = "PUBLIC_KEY";

    std::string serial_str;
    boost::iostreams::back_insert_device<std::string> inserter(serial_str);
    boost::iostreams::stream<boost::iostreams::back_insert_device<std::string>> s(inserter);
    boost::archive::binary_oarchive oa(s);  

    oa << cmd_str << pke;

    s.flush();

    channels_mtx[dstTid].lock();
    channels[dstTid]->asyncSend(std::move(std::unique_ptr<osuCrypto::ByteStream>(new osuCrypto::ByteStream((uint8_t*)serial_str.data(), serial_str.size())))); 
    channels_mtx[dstTid].unlock();    
}

void TaskComm::sendFinish() {
    std::string cmd_str = "FINISH";

    std::string serial_str;
    boost::iostreams::back_insert_device<std::string> inserter(serial_str);
    boost::iostreams::stream<boost::iostreams::back_insert_device<std::string>> s(inserter);
    boost::archive::binary_oarchive oa(s);  

    oa << cmd_str;

    s.flush();

    for (int dstTid=0; dstTid<tileNum_; dstTid++) {
        if (dstTid != tileIndex_) {
            channels_mtx[dstTid].lock();
            channels[dstTid]->asyncSend(std::move(std::unique_ptr<osuCrypto::ByteStream>(new osuCrypto::ByteStream((uint8_t*)serial_str.data(), serial_str.size())))); 
            channels_mtx[dstTid].unlock();  
        }
    } 
}

void TaskComm::recvFinish() {
    for (int srcTid=0; srcTid<tileNum_; srcTid++) {
        if (srcTid != tileIndex_) {
            osuCrypto::ByteStream recv;
            channels[srcTid]->recv(recv);
            std::string serial_str((char*)recv.data(), recv.size());

            // wrap buffer inside a stream and deserialize serial_str into obj
            boost::iostreams::basic_array_source<char> device(serial_str.data(), serial_str.size());
            boost::iostreams::stream<boost::iostreams::basic_array_source<char> > s(device);
            boost::archive::binary_iarchive ia(s);

            std::string cmd_str;
            ia >> cmd_str;

            if (cmd_str != "FINISH") {
                printf("Did not receive expected FINISH!\n");
                exit(-1);
            }
        }
    }
}

std::vector<TaskqHandlerConfig>& TaskComm::getTaskqHandlerConfig() {
    return thc_;
}

TaskqHandlerConfig& TaskComm::getTaskqHandlerConfig(uint32_t tid) {
    return thc_[tid];
}

void TaskComm::doTaskOperandRotation(std::vector<Task>& taskv, uint64_t rotation) {
    if (rotation != 0) {
        int64_t end = (int64_t)(taskv.size() - rotation);
        for (int j=0; j<end; ++j) {
            auto& dst_task = taskv[j];
            auto& src_task = taskv[j + rotation];
            TaskPayload* dst_tp = (TaskPayload*)dst_task.buf;
            TaskPayload* src_tp = (TaskPayload*)src_task.buf;
            dst_tp->copyOperand(*src_tp, 1, 1);
            dst_tp->useMask = true;
            if (dst_task.vertexIndex == src_task.vertexIndex && (!src_task.isDummy)) {
                dst_tp->operandMask = true;
            } else {
                dst_tp->operandMask = false;
            }
        }
    } else {
        for (int j=0; j<taskv.size(); ++j) {
            auto& dst_task = taskv[j];
            TaskPayload* dst_tp = (TaskPayload*)dst_task.buf;
            dst_tp->operandMask = true;
            dst_tp->useMask = false;
        }        
    }    
}

void TaskComm::sendCipherEntryVec(const CipherEntryVec& cev, const uint32_t dstTid) {
    channels_mtx[dstTid].lock();

    std::vector<osuCrypto::Channel*> curChannels(channels);
    auto sendCEV = [this, &curChannels, dstTid](const CipherEntryVec& cur_cev) {
        std::string serial_str;
        boost::iostreams::back_insert_device<std::string> inserter(serial_str);
        boost::iostreams::stream<boost::iostreams::back_insert_device<std::string>> s(inserter);
        boost::archive::binary_oarchive oa(s); 
        std::string send_cmd_str = "CIPHER_ENTRY_VEC";
        oa << send_cmd_str << cur_cev;
        s.flush();
        curChannels[dstTid]->asyncSend(std::move(std::unique_ptr<osuCrypto::ByteStream>(new osuCrypto::ByteStream((uint8_t*)serial_str.data(), serial_str.size()))));        
    };

    auto sendCMD = [this, &curChannels, dstTid](const std::string cur_cmd) {
        std::string serial_str;
        boost::iostreams::back_insert_device<std::string> inserter(serial_str);
        boost::iostreams::stream<boost::iostreams::back_insert_device<std::string>> s(inserter);
        boost::archive::binary_oarchive oa(s);  
        oa << cur_cmd;
        s.flush();
        curChannels[dstTid]->asyncSend(std::move(std::unique_ptr<osuCrypto::ByteStream>(new osuCrypto::ByteStream((uint8_t*)serial_str.data(), serial_str.size()))));        
    };

    std::string begin_cmd_str = "BEGIN_CIPHER_ENTRY_VEC";
    std::string end_cmd_str = "END_CIPHER_ENTRY_VEC";

    sendCMD(begin_cmd_str);

    int64_t step = MAX_SEND_CEV_SIZE;
    int64_t vecSize = cev.ceVecs.size();
    for (int64_t i=0;i<vecSize;i+=step) {
        int64_t endIndex = i+step<vecSize? i+step:vecSize;
        auto first = cev.ceVecs.cbegin() + i;
        auto last = cev.ceVecs.cbegin() + endIndex;
        std::vector<CipherEntry> curVec(first, last);
        CipherEntryVec curCEV;
        curCEV.ceVecs.swap(curVec);

        sendCEV(curCEV);
    }

    sendCMD(end_cmd_str);
 
    channels_mtx[dstTid].unlock(); 
}

void TaskComm::recvCipherEntryVec(CipherEntryVec& cev, const uint32_t srcTid) {
    std::string begin_cmd_str = "BEGIN_CIPHER_ENTRY_VEC";
    std::string end_cmd_str = "END_CIPHER_ENTRY_VEC";
    std::string send_cmd_str = "CIPHER_ENTRY_VEC";
    cev.ceVecs.clear();
    while (true) {
        osuCrypto::ByteStream recv;
        channels[srcTid]->recv(recv);
        std::string serial_str((char*)recv.data(), recv.size());

        // wrap buffer inside a stream and deserialize serial_str into obj
        boost::iostreams::basic_array_source<char> device(serial_str.data(), serial_str.size());
        boost::iostreams::stream<boost::iostreams::basic_array_source<char> > s(device);
        boost::archive::binary_iarchive ia(s);

        std::string cmd_str;
        CipherEntryVec curCEV;
        ia >> cmd_str;

        if (cmd_str == begin_cmd_str) {
            continue;
        } else if (cmd_str == end_cmd_str) {
            break;
        } else if (cmd_str == send_cmd_str) {
            ia >> curCEV;
            cev.ceVecs.insert(cev.ceVecs.end(), curCEV.ceVecs.begin(), curCEV.ceVecs.end());
            continue;
        } else {
            printf("Did not receive expected CIPHER_ENTRY_VEC!\n");
            std::cout<<cmd_str<<std::endl;
            exit(-1);
        }
    }
}

void TaskComm::sendShareVecVec(const ShareVecVec& svv, const uint32_t dstTid) {
    std::string cmd_str = "SHARE_VEC_VEC";

    std::string serial_str;
    boost::iostreams::back_insert_device<std::string> inserter(serial_str);
    boost::iostreams::stream<boost::iostreams::back_insert_device<std::string>> s(inserter);
    boost::archive::binary_oarchive oa(s);  

    oa << cmd_str << svv;

    s.flush();

    channels_mtx[dstTid].lock();
    channels[dstTid]->asyncSend(std::move(std::unique_ptr<osuCrypto::ByteStream>(new osuCrypto::ByteStream((uint8_t*)serial_str.data(), serial_str.size())))); 
    channels_mtx[dstTid].unlock(); 
}

void TaskComm::recvShareVecVec(ShareVecVec& svv, const uint32_t srcTid) {
    osuCrypto::ByteStream recv;
    channels[srcTid]->recv(recv);
    std::string serial_str((char*)recv.data(), recv.size());

    // wrap buffer inside a stream and deserialize serial_str into obj
    boost::iostreams::basic_array_source<char> device(serial_str.data(), serial_str.size());
    boost::iostreams::stream<boost::iostreams::basic_array_source<char> > s(device);
    boost::archive::binary_iarchive ia(s);

    std::string cmd_str;
    ia >> cmd_str >> svv;

    if (cmd_str != "SHARE_VEC_VEC") {
        printf("Did not receive expected SHARE_VEC_VEC!\n");
        exit(-1);
    }
}

void TaskComm::sendObliviousMapperConfig(const uint64_t srcPosNum, const uint64_t dstPosNum, const uint64_t plainNumPerPos, const uint32_t dstTid) {
    std::string cmd_str = "OM_CONFIG";

    std::string serial_str;
    boost::iostreams::back_insert_device<std::string> inserter(serial_str);
    boost::iostreams::stream<boost::iostreams::back_insert_device<std::string>> s(inserter);
    boost::archive::binary_oarchive oa(s);  

    oa << cmd_str << srcPosNum << dstPosNum << plainNumPerPos;

    s.flush();

    channels_mtx[dstTid].lock();
    channels[dstTid]->asyncSend(std::move(std::unique_ptr<osuCrypto::ByteStream>(new osuCrypto::ByteStream((uint8_t*)serial_str.data(), serial_str.size())))); 
    channels_mtx[dstTid].unlock();     
}

void TaskComm::recvObliviousMapperConfig(uint64_t& srcPosNum, uint64_t& dstPosNum, uint64_t& plainNumPerPos, const uint32_t srcTid) {
    osuCrypto::ByteStream recv;
    channels[srcTid]->recv(recv);
    std::string serial_str((char*)recv.data(), recv.size());

    // wrap buffer inside a stream and deserialize serial_str into obj
    boost::iostreams::basic_array_source<char> device(serial_str.data(), serial_str.size());
    boost::iostreams::stream<boost::iostreams::basic_array_source<char> > s(device);
    boost::archive::binary_iarchive ia(s);

    std::string cmd_str;
    ia >> cmd_str >> srcPosNum >> dstPosNum >> plainNumPerPos;

    if (cmd_str != "OM_CONFIG") {
        printf("Did not receive expected OM_CONFIG!\n");
        exit(-1);
    }
}

