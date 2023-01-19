#ifndef SSHE_TASKQ_HANDLER_H_
#define SSHE_TASKQ_HANDLER_H_

#include "task.h"
#include "TaskUtil.h"

#include "Common/Defines.h"
#include "Network/BtIOService.h"

#include "Network/BtEndpoint.h"
#include "Network/Channel.h"

#include "Common/ByteStream.h"

#include "emp-sh2pc/emp-sh2pc.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>
#include <algorithm>
#include <errno.h>
#include <memory>
#include <chrono>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/access.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/back_inserter.hpp>

void print_duration(std::chrono::_V2::system_clock::time_point t1, string tag);
void print_comm(uint64_t c_byte, string tag);

void task_queue_handler(std::queue<Task>* taskq);
void issue_client_recv_threads(std::vector<std::thread>& threads);
void issue_server_recv_threads(std::vector<std::thread>& threads);
// void issue_client_computeTaskq_threads(std::vector<std::thread>& threads);
void split_task_queue_on_destination(std::queue<Task>* taskq, std::vector<std::queue<Task>>& taskqs, std::queue<Task>& plaintext_taskq);
void split_task_queue_on_destination(std::vector<Task>& taskq, std::vector<std::vector<Task>>& taskvs, std::vector<Task>& plaintext_taskv);
void split_task_queue_on_destination(std::queue<Task>* taskq, std::vector<std::vector<Task>>& taskvs, std::vector<Task>& plaintext_taskv);

void single_destination_task_vector_handler(std::vector<Task>& taskv, size_t dstTid);
void plain_task_vector_handler(std::vector<Task>& taskv);
uint64_t get_rotation_upper_bound(std::vector<Task>& taskv);
uint64_t client_task_vector_direct_handler(std::vector<Task>& taskv, size_t dstTid);
void client_direct_recv(std::vector<Task>& taskv, size_t dstTid);
uint64_t server_task_vector_direct_handler(std::vector<Task>& taskv, size_t srcTid);
uint64_t get_next_power_of_2(uint64_t num);

struct TaskqHandlerConfig {
    bool sendOperand;
    bool mergeResult;
    bool sendTaskqDigest;
    uint64_t rotation;
};
typedef struct TaskqHandlerConfig TaskqHandlerConfig;

// Singleton task data communication object
class TaskComm {
public:
    TaskComm(TaskComm const&) = delete;
    TaskComm& operator=(TaskComm const&) = delete;

    /**
     * @brief Get the client instance object. Client instance sends tasks to other parties.
     * 
     * @return TaskComm& 
     */
    static TaskComm& getClientInstance(){
        static TaskComm theInstance;
        return theInstance;
    }

    /**
     * @brief Get the server instance object. Server instance receives tasks from other parties.
     * 
     * @return TaskComm& 
     */
    static TaskComm& getServerInstance(){
        static TaskComm theInstance;
        return theInstance;
    }

    /**
     * Set remote networking channels, public keys and other configurations.
     */
    void setUp(bool isClient);

    /**
     * Clean task queue-related data after computation.
     */
    void clean();

    /**
     * Clean task queue-related data for a specific dstTid after computation.
     */
    void clean(const uint32_t dstTid);

    void closeChannels();

    void tileNumIs(size_t tileNum) {
        tileNum_ = tileNum;
    }

    size_t getTileNum() {
        return tileNum_;
    }

    void tileIndexIs(size_t tileIndex) {
        tileIndex_ = tileIndex;
    }

    void settingIs(std::string setting) {
        setting_ = setting;
    }

    std::string getSetting() {
        return setting_;
    }

    void noPreprocessIs(bool noPreprocess) {
        noPreprocess_ = noPreprocess;
    }

    bool getNoPreprocess() {
        return noPreprocess_;
    }

    void isClusterIs(bool isCluster) {
        isCluster_ = isCluster;
    }

    bool getIsCluster() {
        return isCluster_;
    }

    void isRotationBasedIs(bool isRotationBased) {
        isRotationBased_ = isRotationBased;
    }

    bool getIsRotationBased() {
        return isRotationBased_;
    }

    size_t getTileIndex() {
        return tileIndex_;
    }

    std::vector<std::vector<Task>>& getTaskvs() {
        return taskvs;
    }

    std::vector<std::vector<Task>>& getLocalTaskvs() {
        return local_taskvs;
    }    

    std::vector<Task>& getTaskv(const uint32_t tid) {
        return taskvs[tid];
    }

    std::vector<std::mutex>& getClientComputeTaskqMtx() {
        return client_computeTaskq_mtx;
    }

    std::vector<std::condition_variable>& getClientComputeTaskqCv() {
        return client_computeTaskq_cv;
    }
    
    std::vector<std::atomic<bool>>& getClientWorkdone() {
        return client_workdone;
    }

    std::vector<Semaphore>& getClientComputeTaskqSmp() {
        return client_computeTaskq_smp;
    }

    Semaphore& getRemoteUpdateReadySmp(const uint32_t index) {
        return remote_update_ready_smp[index];
    }

    Semaphore& getLocalUpdateReadySmp(const uint32_t index) {
        return local_update_ready_smp[index];
    }

    uint32_t getMPCBasePort() {
        return mpcBasePort;
    }

    std::string getMPCIP() {
        return mpcIp;
    }

    void sendTaskDigest(Task task, const uint32_t dstTid);
    void sendTaskqDigest(Task task, const uint64_t taskNum, const uint32_t dstTid);
    void sendExpectedTaskNum(const uint64_t taskNum, const uint64_t dstTid);
    void sendDecrypt(CipherEntry& ce, const uint32_t dstTid, const uint32_t taskId, 
                        const uint32_t operandId0);
    void sendDecryptPacked(std::shared_ptr<PackedCipherEntry> pce);
    void sendEndDecrypt(const uint32_t dstTid);

    // Send encrypted share
    void sendEncrypted(CipherEntry& ce, const uint32_t dstTid, const uint32_t taskId, 
                        const uint32_t operandId0);
    void sendEndEncrypted(const uint32_t dstTid);

    int recvEntrance(const uint32_t srcTid);

    void taskDigestHandler(Task task, const uint32_t srcTid);
    void taskqDigestHandler(Task emptyTask, const uint64_t taskNum, const uint32_t srcTid);
    // void writeCipherEntryHandler(CipherEntry& ce, const uint32_t srcTid, const uint32_t dstTid, const uint32_t taskId, const uint32_t operandId0);
    void decryptHandler(CipherEntry& ce, const uint32_t srcTid, const uint32_t dstTid, const uint32_t taskId, const uint32_t operandId0);
    void decryptPackedHandler(const uint32_t srcTid);
    void writeShareHandler(const ShareVec& shareVec, const uint32_t srcTid, const uint32_t taskId, const uint32_t operandId0);
    void mergeEncryptedHandler(CipherEntry& ce, const uint32_t srcTid, const uint32_t taskId, const uint32_t operandId0);
    void sendRouteShare(const ShareVec& shareVec, const uint32_t srcTid, const uint32_t dstTid, const uint32_t taskId, const uint32_t operandId0);
    void sendPlain(const uint32_t dstTid, const uint32_t taskId, const uint32_t operandId0);
    void writePlainHandler(const uint32_t srcTid, const uint32_t taskId, const uint32_t operandId0);
    void computeTaskq(vector<Task>& taskq, uint32_t curTid, uint32_t dstTid, int party, int rotation = 0);
    void computePlaintextTaskq(vector<Task>& taskq);

    void ssAddPairDoubleUint(Task& task, int party);
    void ssAddMixedPairDoubleUint(Task& task, int party);
    void ssMinUintWithParent(Task& task, int party);
    void ssAddDouble(Task& task, int party);
    void ssAddUlong(Task& task, int party);
    void ssAddUint(Task& task, int party);
    void ssAddUintWithReplaceParent(Task& task, int party);
    void ssDivDouble(Task& task, int party);
    void ssSwapCipherEntry(Task& task, int party);

    void sendEncryptedTaskqShare(const uint32_t dstTid);
    void sendEncryptedTaskqShare(const std::vector<Task>& taskv, const uint32_t dstTid);

    void sendObliviousMapperConfig(const uint64_t srcPosNum, const uint64_t dstPosNum, const uint64_t plainNumPerPos, const uint32_t dstTid);
    void recvObliviousMapperConfig(uint64_t& srcPosNum, uint64_t& dstPosNum, uint64_t& plainNumPerPos, const uint32_t srcTid);

    void sendCipherEntryVec(const CipherEntryVec& cev, const uint32_t dstTid);
    void recvCipherEntryVec(CipherEntryVec& cev, const uint32_t srcTid);

    void sendShareVecVec(const ShareVecVec& svv, const uint32_t dstTid);
    void recvShareVecVec(ShareVecVec& svv, const uint32_t srcTid);

    void setExpectedTaskNum(const uint64_t num, const uint32_t dstTid) {
        expectedTaskNum[dstTid] = num;
        if (expectedTaskNum[dstTid] == finishTaskNum[dstTid])
            localFinishqFlag[dstTid] = true;
    }

    void sendRotationUB(const uint64_t rotation_ub, const uint32_t dstTid);

    void setRemoteFinishqFlag(const uint32_t dstTid);
    bool getRemoteFinishqFlag(const uint32_t dstTid);
    void setLocalFinishqFlag(const uint32_t srcTid);
    bool getLocalFinishqFlag(const uint32_t srcTid);

    std::vector<std::thread>& getComputeRemoteTaskqThreads();

    void sendPublicKey(const uint32_t dstTid, const PublicKeyEntry& pke);
    void sendFinish();
    void recvFinish();
    void sendEndRotation(const uint32_t dstTid);
    std::vector<TaskqHandlerConfig>& getTaskqHandlerConfig();
    TaskqHandlerConfig& getTaskqHandlerConfig(uint32_t tid);
    void doTaskOperandRotation(std::vector<Task>& taskv, uint64_t rotation);

private:
    // Remote networking channel
    std::vector<osuCrypto::Channel*> channels;
    std::vector<std::mutex> channels_mtx;
    std::vector<std::mutex> client_computeTaskq_mtx;
    std::vector<std::condition_variable> client_computeTaskq_cv;
    std::vector<std::atomic<bool>> client_workdone;
    std::vector<Semaphore> client_computeTaskq_smp;
    std::vector<Semaphore> remote_update_ready_smp;
    std::vector<Semaphore> local_update_ready_smp;
    uint32_t basePort;
    std::string ip;
    uint32_t mpcBasePort;
    std::string mpcIp;
    std::list<osuCrypto::BtEndpoint> endpoints;
    osuCrypto::BtIOService* ioService;
    size_t tileNum_;
    size_t tileIndex_;
    std::vector<std::vector<Task>> taskvs; // Remote (received) task
    std::vector<std::vector<Task>> local_taskvs; // Local task
    std::vector<uint64_t> expectedTaskNum;
    std::vector<std::atomic<uint64_t>> readyTaskNum;
    std::vector<std::atomic<uint64_t>> finishTaskNum;
    CryptoUtil& cryptoUtil;
    std::vector<bool> remoteFinishqFlag;
    std::vector<bool> localFinishqFlag;
    std::vector<std::thread> computeRemoteTaskqThreads;
    // std::vector<std::thread> decryptThreads;
    std::vector<std::vector<Task>> decrypt_taskq; // for parallization
    std::vector<std::vector<Task>> merge_encrypted_taskq; // for parallization
    std::vector<std::vector<std::shared_ptr<PackedCipherEntry>>> packed_ce_q; // packed cipher entry queue
    bool isClient_;
    std::vector<TaskqHandlerConfig> thc_;
    std::string setting_;
    bool noPreprocess_;
    bool isCluster_;
    bool isRotationBased_;

    TaskComm(): cryptoUtil(CryptoUtil::getInstance()) {}
};

#endif