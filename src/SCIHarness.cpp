#include "SCIHarness.h"

namespace sci{

#define print_fixed(vec)                                                         \
  {                                                                            \
    auto tmp_pub = fix_op->output(PUBLIC, vec).subset(I, I + J);                  \
    std::cout << #vec << "_pub: " << tmp_pub << std::endl;                               \
  }

template<typename T>
std::vector<T> slice(std::vector<T> const &v, int m, int n) {
   auto first = v.begin() + m;
   auto last = v.begin() + n;
   std::vector<T> vector(first, last);
   return vector;
}

uint64_t computeTaskq(std::vector<Task>& taskq, uint32_t curTid, uint32_t dstTid, uint32_t mpcBasePort, uint32_t tileNum, uint32_t tileIndex, int party, int rotation) {
    auto t1 = std::chrono::high_resolution_clock::now();

    TaskComm& taskComm = TaskComm::getClientInstance();
    bool isCluster = taskComm.getIsCluster();

    printf("Begin SCI compute taskq length %d %d\n", taskq.size(), party);
    printf("rotation %d\n", rotation);
    int64_t cur_taskv_size = taskq.size() - rotation;
    std::vector<Task> taskv = slice(taskq, 0, cur_taskv_size);

    // print("cur_taskv_size %d\n")
    if (cur_taskv_size <= 0) {
        printf("End compute taskq length %d %d\n", taskq.size(), party);
        return 0;
    }
    
    IOPack *iopack = nullptr;
    OTPack *otpack = nullptr;
    FixOp *fix_op = nullptr;
    FPOp *fp_op = nullptr;
    int port = party==emp::ALICE? (mpcBasePort + dstTid * tileNum + tileIndex) : (mpcBasePort + tileIndex * tileNum + dstTid);

    string address = "";
    if (!isCluster) {
        address = "127.0.0.1";
    } else {
        address = "10.0.0.";
        if (party==emp::ALICE) address += std::to_string(tileIndex + 1);
        else address += std::to_string(dstTid + 1);
    }

    iopack = new IOPack(party, port, address);
    otpack = new OTPack(iopack, party);
    fix_op = new FixOp(party, iopack, otpack);
    fp_op = new FPOp(party, iopack, otpack);

    uint64_t comm_start = iopack->get_comm();   

    Task& task_0 = taskv[0];
    switch (task_0.type) {
        case ADD_PAIR_DOUBLE_UINT: {
            ssAddPairDoubleUint(taskv, fix_op, fp_op, party);
        } break;
        case ADD_MIXED_PAIR_DOUBLE_UINT: {
            ssAddMixedPairDoubleUint(taskv, fix_op, fp_op, party);
        } break;
        case MIN_UINT_WITH_PARENT: {
            ssMinUintWithParent(taskv, fix_op, fp_op, party);
        } break;
        case ADD_DOUBLE: {
            ssAddDouble(taskv, fix_op, fp_op, party);
        } break;
        case ADD_ULONG: {
            ssAddUlong(taskv, fix_op, fp_op, party);
        } break;
        case ADD_UINT: {
            ssAddUint(taskv, fix_op, fp_op, party);
        } break;
        case ADD_UINT_WITH_REPLACE_PARENT: {
            ssAddUintWithReplaceParent(taskv, fix_op, fp_op, party);
        } break;
        case UINT_REPLACE_PARENT: {
            ssUintReplaceParent(taskv, fix_op, fp_op, party);
        } break;
        case DIV_DOUBLE: {
            ssDivDouble(taskv, fix_op, fp_op, party);
        } break;
        default: {
        } break;
    }
    
    // tp->isOperandEncrypted[0] = true;
    // if (rotation != 0) {
    //     tp->unifyOperand();
    // }
    uint64_t comm_end = iopack->get_comm();

    delete iopack;
    delete otpack;
    delete fix_op;
    delete fp_op;

    printf("End compute taskq length %d %d\n", taskq.size(), party);

    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "::computetaskq took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() / 1000
              << " seconds\n";
    return (comm_end - comm_start);
}

void ssAddPairDoubleUint(std::vector<Task>& taskv, FixOp* fix_op, FPOp* fp_op, int party) {

}

void ssAddMixedPairDoubleUint(std::vector<Task>& taskv, FixOp* fix_op, FPOp* fp_op, int party) {

}

void ssMinUintWithParent(std::vector<Task>& taskv, FixOp* fix_op, FPOp* fp_op, int party) {
    auto task_content0 = (MinWithParent<uint32_t, uint64_t>*)taskv[0].buf;
    uint64_t taskNum = taskv.size();

	struct {
		FixArray val;
		FixArray parent;
	} sci_operands[3];

    struct {
        uint64_t* val;
        uint64_t* parent;
    } operands_share[2];

    uint8_t* operandMask = new uint8_t[taskNum];

    for (int j=0; j<taskNum; ++j) {
        auto task_content = (AddWithReplaceParent<uint32_t, uint64_t>*)taskv[j].buf;
        if (task_content->operandMask && (!taskv[j].isDummy)) operandMask[j] = 1;
        else operandMask[j] = 0;
    }
    BoolArray sci_operandMask = fix_op->bool_op->input(ALICE, taskNum, operandMask);

    for (int i=0; i<2; ++i) {
        operands_share[i].val = new uint64_t[taskNum];
        operands_share[i].parent = new uint64_t[taskNum];

        if (!task_content0->isOperandEncrypted[i]) {
            for (int j=0; j<taskNum; ++j) {
                auto task_content = (AddWithReplaceParent<uint32_t, uint64_t>*)taskv[j].buf;
                operands_share[i].val[j] = (uint64_t)task_content->operands[i].plain_val;
                operands_share[i].parent[j] = (uint64_t)task_content->operands[i].plain_parent;
            }
            sci_operands[i].val = fix_op->input(ALICE, taskNum, operands_share[i].val, false, 64);
            sci_operands[i].parent = fix_op->input(ALICE, taskNum, operands_share[i].parent, false, 64);
        } else {
            for (int j=0; j<taskNum; ++j) {
                auto task_content = (AddWithReplaceParent<uint32_t, uint64_t>*)taskv[j].buf;
                operands_share[i].val[j] = (uint64_t)task_content->operands[i].enc.share[0];
                operands_share[i].parent[j] = (uint64_t)task_content->operands[i].enc.share[1];
            }
            sci_operands[i].val = fix_op->input(party, taskNum, operands_share[i].val, false, 64);
            sci_operands[i].parent = fix_op->input(party, taskNum, operands_share[i].parent, false, 64);       
        }
    }

    // Compute
    BoolArray selector;
    selector = fix_op->GT(sci_operands[0].val, sci_operands[1].val);
    sci_operands[2].val = fix_op->if_else(selector, sci_operands[1].val, sci_operands[0].val);
    sci_operands[2].parent = fix_op->if_else(selector, sci_operands[1].parent, sci_operands[0].parent);

    // Mask
    if (task_content0->useMask) {
        sci_operands[0].val = fix_op->if_else(sci_operandMask, sci_operands[2].val, sci_operands[0].val);
        sci_operands[0].parent = fix_op->if_else(sci_operandMask, sci_operands[2].parent, sci_operands[0].parent);
    } else {
        sci_operands[0].val = sci_operands[2].val;
        sci_operands[0].parent = sci_operands[2].parent;
    }

    // Write back results
    for (int j=0; j<taskNum; ++j) {
        auto task_content = (AddWithReplaceParent<uint32_t, uint64_t>*)taskv[j].buf;
        task_content->operands[0].enc.share[0] = (int64_t)sci_operands[0].val.data[j];
        task_content->operands[0].enc.share[1] = (int64_t)sci_operands[0].parent.data[j];
        task_content->unifyOperand();
    }

    for (int i=0; i<2; ++i) {
        delete operands_share[i].val;
        delete operands_share[i].parent;
    }
    delete operandMask;
}

void ssAddDouble(std::vector<Task>& taskv, FixOp* fix_op, FPOp* fp_op, int party) {
    auto task_content0 = (Add<double>*)taskv[0].buf;
    uint64_t taskNum = taskv.size();

    FixArray sci_operands[3];

    uint64_t* operands_share[2];

    uint8_t* operandMask = new uint8_t[taskNum];

    for (int j=0; j<taskNum; ++j) {
        auto task_content = (Add<double>*)taskv[j].buf;
        if (task_content->operandMask && (!taskv[j].isDummy)) operandMask[j] = 1;
        else operandMask[j] = 0;
    }
    BoolArray sci_operandMask = fix_op->bool_op->input(ALICE, taskNum, operandMask);

    for (int i=0; i<2; ++i) {
        operands_share[i] = new uint64_t[taskNum];

        if (!task_content0->isOperandEncrypted[i]) {
            for (int j=0; j<taskNum; ++j) {
                auto task_content = (Add<double>*)taskv[j].buf;
                operands_share[i][j] = (uint64_t)(task_content->operands[i].plain * (1<<16));
            }
            sci_operands[i] = fix_op->input(ALICE, taskNum, operands_share[i], false, 64, 16);
        } else {
            for (int j=0; j<taskNum; ++j) {
                auto task_content = (Add<double>*)taskv[j].buf;
                operands_share[i][j] = (uint64_t)task_content->operands[i].enc.share[0];
            }
            sci_operands[i] = fix_op->input(party, taskNum, operands_share[i], false, 64, 16);
        }
    }

    // Compute
    sci_operands[2] = fix_op->add(sci_operands[0], sci_operands[1]);

    // Mask
    if (task_content0->useMask) {
        sci_operands[0] = fix_op->if_else(sci_operandMask, sci_operands[2], sci_operands[0]);
    } else {
        sci_operands[0] = sci_operands[2];
    }

    // print_fixed(sci_operands[0]);

    // Write back results
    for (int j=0; j<taskNum; ++j) {
        auto task_content = (Add<double>*)taskv[j].buf;
        task_content->operands[0].enc.share[0] = (int64_t)sci_operands[0].data[j];
        // if (party == ALICE) std::cout<<"Gather share: "<<task_content->operands[0].enc.share[0]<<std::endl;
        task_content->unifyOperand();
    }

    for (int i=0; i<2; ++i) {
        delete operands_share[i];
    }
    delete operandMask;
}

void ssAddUlong(std::vector<Task>& taskv, FixOp* fix_op, FPOp* fp_op, int party) {

}

void ssAddUint(std::vector<Task>& taskv, FixOp* fix_op, FPOp* fp_op, int party) {

}

void ssAddUintWithReplaceParent(std::vector<Task>& taskv, FixOp* fix_op, FPOp* fp_op, int party) {
    auto task_content0 = (AddWithReplaceParent<uint32_t, uint64_t>*)taskv[0].buf;
    uint64_t taskNum = taskv.size();

	struct {
		FixArray val;
		FixArray parent;
	} sci_operands[3];

    struct {
        uint64_t* val;
        uint64_t* parent;
    } operands_share[2];

    uint8_t* operandMask = new uint8_t[taskNum];

    for (int j=0; j<taskNum; ++j) {
        auto task_content = (AddWithReplaceParent<uint32_t, uint64_t>*)taskv[j].buf;
        if (task_content->operandMask && (!taskv[j].isDummy)) operandMask[j] = 1;
        else operandMask[j] = 0;
    }
    BoolArray sci_operandMask = fix_op->bool_op->input(ALICE, taskNum, operandMask);

    for (int i=0; i<2; ++i) {
        operands_share[i].val = new uint64_t[taskNum];
        operands_share[i].parent = new uint64_t[taskNum];

        if (!task_content0->isOperandEncrypted[i]) {
            for (int j=0; j<taskNum; ++j) {
                auto task_content = (AddWithReplaceParent<uint32_t, uint64_t>*)taskv[j].buf;
                operands_share[i].val[j] = (uint64_t)task_content->operands[i].plain_val;
                operands_share[i].parent[j] = (uint64_t)task_content->operands[i].plain_parent;
            }
            sci_operands[i].val = fix_op->input(ALICE, taskNum, operands_share[i].val, false, 64);
            sci_operands[i].parent = fix_op->input(ALICE, taskNum, operands_share[i].parent, false, 64);
        } else {
            for (int j=0; j<taskNum; ++j) {
                auto task_content = (AddWithReplaceParent<uint32_t, uint64_t>*)taskv[j].buf;
                operands_share[i].val[j] = (uint64_t)task_content->operands[i].enc.share[0];
                operands_share[i].parent[j] = (uint64_t)task_content->operands[i].enc.share[1];
            }
            sci_operands[i].val = fix_op->input(party, taskNum, operands_share[i].val, false, 64);
            sci_operands[i].parent = fix_op->input(party, taskNum, operands_share[i].parent, false, 64);       
        }
    }

    // Compute
    sci_operands[2].val = fix_op->add(sci_operands[0].val, sci_operands[1].val);
    sci_operands[2].parent = sci_operands[1].parent;

    // Mask
    if (task_content0->useMask) {
        sci_operands[0].val = fix_op->if_else(sci_operandMask, sci_operands[2].val, sci_operands[0].val);
        sci_operands[0].parent = fix_op->if_else(sci_operandMask, sci_operands[2].parent, sci_operands[0].parent);
    } else {
        sci_operands[0].val = sci_operands[2].val;
        sci_operands[0].parent = sci_operands[2].parent;
    }

    // Write back results
    for (int j=0; j<taskNum; ++j) {
        auto task_content = (AddWithReplaceParent<uint32_t, uint64_t>*)taskv[j].buf;
        task_content->operands[0].enc.share[0] = (int64_t)sci_operands[0].val.data[j];
        task_content->operands[0].enc.share[1] = (int64_t)sci_operands[0].parent.data[j];
        task_content->unifyOperand();
    }

    for (int i=0; i<2; ++i) {
        delete operands_share[i].val;
        delete operands_share[i].parent;
    }
    delete operandMask;
}

void ssUintReplaceParent(std::vector<Task>& taskv, FixOp* fix_op, FPOp* fp_op, int party) {
    auto task_content0 = (ReplaceParent<uint32_t, uint64_t>*)taskv[0].buf;
    uint64_t taskNum = taskv.size();

	struct {
		FixArray val;
		FixArray parent;
	} sci_operands[3];

    struct {
        uint64_t* val;
        uint64_t* parent;
    } operands_share[2];

    uint8_t* operandMask = new uint8_t[taskNum];

    for (int j=0; j<taskNum; ++j) {
        auto task_content = (ReplaceParent<uint32_t, uint64_t>*)taskv[j].buf;
        if (task_content->operandMask && (!taskv[j].isDummy)) operandMask[j] = 1;
        else operandMask[j] = 0;
    }
    BoolArray sci_operandMask = fix_op->bool_op->input(ALICE, taskNum, operandMask);

    for (int i=0; i<2; ++i) {
        operands_share[i].val = new uint64_t[taskNum];
        operands_share[i].parent = new uint64_t[taskNum];

        if (!task_content0->isOperandEncrypted[i]) {
            for (int j=0; j<taskNum; ++j) {
                auto task_content = (ReplaceParent<uint32_t, uint64_t>*)taskv[j].buf;
                operands_share[i].val[j] = (uint64_t)task_content->operands[i].plain_val;
                operands_share[i].parent[j] = (uint64_t)task_content->operands[i].plain_parent;
            }
            sci_operands[i].val = fix_op->input(ALICE, taskNum, operands_share[i].val, false, 64);
            sci_operands[i].parent = fix_op->input(ALICE, taskNum, operands_share[i].parent, false, 64);
        } else {
            for (int j=0; j<taskNum; ++j) {
                auto task_content = (ReplaceParent<uint32_t, uint64_t>*)taskv[j].buf;
                operands_share[i].val[j] = (uint64_t)task_content->operands[i].enc.share[0];
                operands_share[i].parent[j] = (uint64_t)task_content->operands[i].enc.share[1];
            }
            sci_operands[i].val = fix_op->input(party, taskNum, operands_share[i].val, false, 64);
            sci_operands[i].parent = fix_op->input(party, taskNum, operands_share[i].parent, false, 64);       
        }
    }

    // Compute
    sci_operands[2].val = sci_operands[0].val;
    sci_operands[2].parent = sci_operands[1].parent;

    // Mask
    if (task_content0->useMask) {
        sci_operands[0].val = fix_op->if_else(sci_operandMask, sci_operands[2].val, sci_operands[0].val);
        sci_operands[0].parent = fix_op->if_else(sci_operandMask, sci_operands[2].parent, sci_operands[0].parent);
    } else {
        sci_operands[0].val = sci_operands[2].val;
        sci_operands[0].parent = sci_operands[2].parent;
    }

    // Write back results
    for (int j=0; j<taskNum; ++j) {
        auto task_content = (ReplaceParent<uint32_t, uint64_t>*)taskv[j].buf;
        task_content->operands[0].enc.share[0] = (int64_t)sci_operands[0].val.data[j];
        task_content->operands[0].enc.share[1] = (int64_t)sci_operands[0].parent.data[j];
        task_content->unifyOperand();
    }

    for (int i=0; i<2; ++i) {
        delete operands_share[i].val;
        delete operands_share[i].parent;
    }
    delete operandMask;    
}

void ssDivDouble(std::vector<Task>& taskv, FixOp* fix_op, FPOp* fp_op, int party) {
    auto task_content0 = (Div<double>*)taskv[0].buf;
    uint64_t taskNum = taskv.size();

    FixArray sci_operands[3];

    uint64_t* operands_share[2];

    uint8_t* operandMask = new uint8_t[taskNum];

    for (int j=0; j<taskNum; ++j) {
        auto task_content = (Div<double>*)taskv[j].buf;
        if (task_content->operandMask && (!taskv[j].isDummy)) operandMask[j] = 1;
        else operandMask[j] = 0;
    }
    BoolArray sci_operandMask = fix_op->bool_op->input(ALICE, taskNum, operandMask);

    for (int i=0; i<2; ++i) {
        operands_share[i] = new uint64_t[taskNum];

        if (!task_content0->isOperandEncrypted[i]) {
            for (int j=0; j<taskNum; ++j) {
                auto task_content = (Div<double>*)taskv[j].buf;
                operands_share[i][j] = (uint64_t)(task_content->operands[i].plain * (1<<16));
            }
            if (i == 0) sci_operands[i] = fix_op->input(ALICE, taskNum, operands_share[i], true, 64, 16); // Risky
            else sci_operands[i] = fix_op->input(ALICE, taskNum, operands_share[i], false, 64, 16); // Risky
            sci_operands[i] = fix_op->reduce(sci_operands[i], 32);
        } else {
            for (int j=0; j<taskNum; ++j) {
                auto task_content = (Div<double>*)taskv[j].buf;
                operands_share[i][j] = (uint64_t)task_content->operands[i].enc.share[0];
            }
            if (i == 0) sci_operands[i] = fix_op->input(party, taskNum, operands_share[i], true, 64, 16); // Risky
            else sci_operands[i] = fix_op->input(party, taskNum, operands_share[i], false, 64, 16);
            sci_operands[i] = fix_op->reduce(sci_operands[i], 32);
        }
    }

    // Compute
    sci_operands[2] = fix_op->div(sci_operands[0], sci_operands[1], 32, 16);

    // Extend
    sci_operands[0] = fix_op->extend(sci_operands[0], 64);
    sci_operands[2] = fix_op->extend(sci_operands[2], 64);

    // Mask
    if (task_content0->useMask) {
        sci_operands[0] = fix_op->if_else(sci_operandMask, sci_operands[2], sci_operands[0]);
    } else {
        sci_operands[0] = sci_operands[2];
    }

    // Write back results
    for (int j=0; j<taskNum; ++j) {
        auto task_content = (Div<double>*)taskv[j].buf;
        task_content->operands[0].enc.share[0] = (int64_t)sci_operands[0].data[j];
        task_content->unifyOperand();
    }

    for (int i=0; i<2; ++i) {
        delete operands_share[i];
    }
    delete operandMask;
}

} // namespace sci