#ifndef OBLIVIOUS_MAPPER_H_
#define OBLIVIOUS_MAPPER_H_

#include "task.h"
#include "TaskqHandler.h"
#include "TaskUtil.h"
#include <string>

void client_oblivious_mapper_preprocess(std::vector<uint64_t>& srcPos, std::vector<uint64_t>& dstPos, 
                                    uint64_t plainNumPerPos, uint64_t iter, uint64_t preprocessId, const uint32_t dstTid, bool allowMissing=false);
uint64_t client_batch_oblivious_mapper_preprocess(std::vector<uint64_t>& srcPos, std::vector<uint64_t>& dstPos, 
                                    uint64_t plainNumPerPos, uint64_t iter, uint64_t preprocessId, const uint32_t dstTid, bool allowMissing=false);

void server_oblivious_mapper_preprocess(uint64_t iter, uint64_t preprocessId, const uint32_t dstTid);
uint64_t server_batch_oblivious_mapper_preprocess(uint64_t iter, uint64_t preprocessId, const uint32_t dstTid);

void client_oblivious_mapper_online(std::vector<uint64_t>& srcPos, std::vector<uint64_t>& dstPos,
                                    const ShareVecVec& srcSvv, ShareVecVec& dstSvv,
                                    uint64_t plainNumPerPos, uint64_t iter, uint64_t preprocessId, const uint32_t dstTid, bool allowMissing=false);

void server_oblivious_mapper_online(const ShareVecVec& srcSvv, ShareVecVec& dstSvv,
                                    uint64_t iter, uint64_t preprocessId, const uint32_t dstTid);

void writeShareVecVectoFile(ShareVecVec& svv, std::string fileName);
void loadShareVecVecFromFile(ShareVecVec& svv, std::string fileName);

#endif