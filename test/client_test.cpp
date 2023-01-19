#include "TaskqHandler.h"
#include "TaskUtil.h"
#include "ObliviousMapper.h"

#include <stdio.h>
#include <errno.h>   // for errno
#include <limits.h>  // for INT_MAX, INT_MIN
#include <stdlib.h>  // for strtol
#include <thread>
#include <vector>
#include <iostream>
#include <thread>
#include <random>
#include <chrono>
#include <algorithm>

int main(int argc, char* argv[]) {
    int threadCount = (int)strtol(argv[1], NULL, 10);
    int tileIndex = (int)strtol(argv[2], NULL, 10);

    TaskComm& clientTaskComm = TaskComm::getClientInstance();
    clientTaskComm.tileNumIs(threadCount);
    clientTaskComm.tileIndexIs(tileIndex);

    printf("HERE1\n");

    TaskComm& serverTaskComm = TaskComm::getServerInstance();
    serverTaskComm.tileNumIs(threadCount);
    serverTaskComm.tileIndexIs(tileIndex);

    printf("HERE2\n");

    CryptoUtil& cryptoUtil = CryptoUtil::getInstance();
    cryptoUtil.tileIndexIs(tileIndex);
    cryptoUtil.setUpPaillierCipher();

    printf("HERE3\n");

    std::thread clientSetupThread([&clientTaskComm](){
        clientTaskComm.setUp(true);
    });

    std::thread serverSetupThread([&serverTaskComm](){
        serverTaskComm.setUp(false);
    });

    clientSetupThread.join();
    serverSetupThread.join();

    printf("HERE4\n");

    std::vector<uint64_t> srcPos;
    std::vector<uint64_t> copySrcPos;
    std::vector<uint64_t> dstPos;

    uint64_t srcPosNum = 200000;
    uint64_t dstPosNum = 100000;
    uint64_t plainNumPerPos = 2;

    for (int i=0; i<srcPosNum; ++i) {
        srcPos.push_back(i+1);
        copySrcPos.push_back(i+1);
    }

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    shuffle(copySrcPos.begin(), copySrcPos.end(), std::default_random_engine(seed));

    for (int i=0; i<dstPosNum; ++i) {
        dstPos.push_back(copySrcPos[i]);
    }

    auto t1 = std::chrono::high_resolution_clock::now();

    // client_oblivious_mapper_preprocess(srcPos, dstPos, plainNumPerPos, 0, 0, 1);

    uint64_t batchSize = client_batch_oblivious_mapper_preprocess(srcPos, dstPos, plainNumPerPos, 0, 0, 1);
    std::cout<< "batch size = "<< batchSize << std::endl;

    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "::Client preprocessing took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() / 1000
              << " seconds\n";

    ShareVecVec srcSvv;
    srcSvv.shareVecs.resize(srcPosNum);
    for (int i=0; i<srcPosNum; ++i) {
        for (int j=0; j<plainNumPerPos; ++j) {
            srcSvv.shareVecs[i].shares.push_back(0);
        }
    }

    ShareVecVec dstSvv;

    t1 = std::chrono::high_resolution_clock::now();

    client_oblivious_mapper_online(srcPos, dstPos, srcSvv, dstSvv, plainNumPerPos, 5, 0, 1);

    t2 = std::chrono::high_resolution_clock::now();
    std::cout << "::Client online took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() / 1000
              << " seconds\n";

    uint64_t printNum = 10;
    for (int i=0; i<printNum; ++i) {
        std::cout<<dstPos[i]<<" "<<dstSvv.shareVecs[i].shares[0]<<" "<<dstSvv.shareVecs[i].shares[1]<<std::endl;
    }

    return 0;
}