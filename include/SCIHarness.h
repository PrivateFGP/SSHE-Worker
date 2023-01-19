#ifndef SCI_HARNESS_H_
#define SCI_HARNESS_H_

#include "task.h"
#include "TaskUtil.h"
#include "TaskqHandler.h"

#include "FloatingPoint/fixed-point.h"
#include "FloatingPoint/floating-point.h"

namespace sci{

uint64_t computeTaskq(vector<Task>& taskq, uint32_t curTid, uint32_t dstTid, uint32_t mpcBasePort, uint32_t tileNum, uint32_t tileIndex, int party, int rotation = 0);

void ssAddPairDoubleUint(std::vector<Task>& taskv, FixOp* fix_op, FPOp* fp_op, int party);
void ssAddMixedPairDoubleUint(std::vector<Task>& taskv, FixOp* fix_op, FPOp* fp_op, int party);
void ssMinUintWithParent(std::vector<Task>& taskv, FixOp* fix_op, FPOp* fp_op, int party);
void ssAddDouble(std::vector<Task>& taskv, FixOp* fix_op, FPOp* fp_op, int party);
void ssAddUlong(std::vector<Task>& taskv, FixOp* fix_op, FPOp* fp_op, int party);
void ssAddUint(std::vector<Task>& taskv, FixOp* fix_op, FPOp* fp_op, int party);
void ssAddUintWithReplaceParent(std::vector<Task>& taskv, FixOp* fix_op, FPOp* fp_op, int party);
void ssUintReplaceParent(std::vector<Task>& taskv, FixOp* fix_op, FPOp* fp_op, int party);
void ssDivDouble(std::vector<Task>& taskv, FixOp* fix_op, FPOp* fp_op, int party);

} // namespace sci


#endif