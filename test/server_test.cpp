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

    auto t1 = std::chrono::high_resolution_clock::now();

    // server_oblivious_mapper_preprocess(0, 0, 0);
    uint64_t batchSize = server_batch_oblivious_mapper_preprocess(0, 0, 0);
    std::cout<< "batch size = "<< batchSize << std::endl;

    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "::Server preprocessing took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() / 1000
              << " seconds\n";

    uint64_t srcPosNum = 200000;
    uint64_t dstPosNum = 0;
    uint64_t plainNumPerPos = 2;

    ShareVecVec srcSvv;
    srcSvv.shareVecs.resize(srcPosNum);
    for (int i=0; i<srcPosNum; ++i) {
        for (int j=0; j<plainNumPerPos; ++j) {
            srcSvv.shareVecs[i].shares.push_back(2*i);
        }
    }

    ShareVecVec dstSvv;

    t1 = std::chrono::high_resolution_clock::now();

    server_oblivious_mapper_online(srcSvv, dstSvv, 5, 0, 0);

    t2 = std::chrono::high_resolution_clock::now();
    std::cout << "::Server online took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() / 1000
              << " seconds\n";

    dstPosNum = dstSvv.shareVecs.size();

    uint64_t printNum = 10;
    for (int i=0; i<printNum; ++i) {
        std::cout<<dstSvv.shareVecs[i].shares[0]<<" "<<dstSvv.shareVecs[i].shares[1]<<std::endl;
    }

    return 0;
}