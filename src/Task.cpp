#include "task.h"
#include "TaskUtil.h"

CryptoUtil& TaskPayload::cryptoUtil = CryptoUtil::getInstance();

CipherEntry SwapCipherEntry::splitShareFromEncryptedOperand(int operandId0) {
    int64_t share = cryptoUtil.splitRandomShareFromCipherEntry(operands[operandId0].enc);
    operands[operandId0].enc.share[0] = share;
    operands[operandId0].enc.isShare = true;
    return operands[operandId0].enc;
}

CipherEntry SwapCipherEntry::encryptShare(int operandId0, const uint32_t tid) {
    CipherEntry ce = cryptoUtil.encryptToCipherEntry(operands[operandId0].enc.share[0], tid);
    return ce;
}

void SwapCipherEntry::mergeEncryptedShare(CipherEntry& ce, int operandId0) {
    cryptoUtil.mergeShareIntoCipherEntry(ce, operands[operandId0].enc.share[0], ce.tid);
    operands[operandId0].enc = ce;
    operands[operandId0].enc.isShare = false;
    isOperandEncrypted[operandId0] = true;
    return;
}

void SwapCipherEntry::writeShareToOperand(int64_t share, int operandId0, int operandId1) {
    operands[operandId0].enc.share[0] = share;
    operands[operandId0].enc.isShare = true;
    return;
}

ShareVec SwapCipherEntry::getOperandShare(int operandId0) {
    ShareVec sv;
    sv.shares.push_back(operands[operandId0].enc.share[0]);
    return sv;
}

int SwapCipherEntry::getEncTid(int operandId0) {
    return operands[operandId0].enc.tid;
}

void SwapCipherEntry::writeEncryptedOperand(CipherEntry& ce, int operandId0) {
    operands[operandId0].enc = ce;
    operands[operandId0].enc.isShare = false;
    isOperandEncrypted[operandId0] = true;
}

CipherEntry* SwapCipherEntry::getCipherEntryPtr(int operandId0) {
    return &operands[operandId0].enc;
}

void SwapCipherEntry::unifyOperand() {
    operands[1] = operands[0];
    isOperandEncrypted[1] = isOperandEncrypted[0];
}

void SwapCipherEntry::copyOperand(TaskPayload& srcTP, int srcOperandId0, int dstOperandId0) {
    operands[dstOperandId0] = static_cast<SwapCipherEntry&>(srcTP).operands[srcOperandId0];
}

void SwapCipherEntry::setCipherEntryPlainNum(int operandId0, int num) {
    operands[operandId0].enc.plainNum = num;
}