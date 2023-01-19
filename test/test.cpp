#include "TaskqHandler.h"
#include "TaskUtil.h"

#include <stdio.h>
#include <errno.h>   // for errno
#include <limits.h>  // for INT_MAX, INT_MIN
#include <stdlib.h>  // for strtol
#include <thread>
#include <vector>
#include <iostream>
#include <thread>

int main(int argc, char* argv[]) {
    // int threadCount = (int)strtol(argv[1], NULL, 10);
    // int tileIndex = (int)strtol(argv[2], NULL, 10);

    // TaskComm& clientTaskComm = TaskComm::getClientInstance();
    // clientTaskComm.tileNumIs(threadCount);
    // clientTaskComm.tileIndexIs(tileIndex);
    // clientTaskComm.setUp(true);

    // printf("HERE1\n");

    // std::vector<std::thread> recv_threads;
	// issue_server_recv_threads(recv_threads);

    // TaskComm& serverTaskComm = TaskComm::getServerInstance();
    // serverTaskComm.tileNumIs(threadCount);
    // serverTaskComm.tileIndexIs(tileIndex);
    // serverTaskComm.setUp(false);

    // printf("HERE2\n");

    CryptoUtil& cryptoUtil = CryptoUtil::getInstance();
    cryptoUtil.tileIndexIs(0);
    cryptoUtil.setUpPaillierCipher();

    // printf("HERE3\n");

    std::vector<int64_t> plain = {47483647, 10000};
    std::vector<int64_t> shares;
    CipherEntry ce = cryptoUtil.encryptToCipherEntry(plain, 0);
    cryptoUtil.splitRandomShareFromCipherEntry(ce, 2, shares);
    printf("%lld, %lld\n", shares[0], shares[1]);

    // cryptoUtil.mergeShareIntoCipherEntry(ce, shares, 0);
    
    std::vector<int64_t> remained_shares;
    cryptoUtil.decryptCipherEntryToInt64Vec(ce, remained_shares);
    // printf("%lld, %lld\n", remained_shares[0], remained_shares[1]);

    printf("%lld, %lld\n", (int64_t)(shares[0] + remained_shares[0]), (int64_t)(shares[1] + remained_shares[1]));
    // ophelib::Integer packed_share = cryptoUtil.shareVecToOPHEInteger(shares);
    // ophelib::Integer packed_remained = cryptoUtil.shareVecToOPHEInteger(remained_shares);
    // ophelib::Integer packed_plain = packed_share + packed_remained;
    // std::cout<< packed_share.to_string_() << std::endl;
    // std::cout<< packed_remained.to_string_() << std::endl;
    // std::vector<int64_t> merged_plain;
    // cryptoUtil.OPHEIntegerToInt64Vec(merged_plain, packed_plain, 2);
    // std::cout<< packed_plain.to_string_() << std::endl;
    // printf("%lld, %lld\n", merged_plain[0], merged_plain[1]);

    // std::cout<<(ophelib::Integer(-1)%ophelib::Integer(100)).to_string_()<<std::endl;

    // std::cout<<ophelib::Integer(-1)<<std::endl;


    std::vector<std::thread> threads;
	for (int i=0; i<5; i++) {
        threads.emplace_back([i] () {
            printf("halo %d\n", i);
        });
	}
    for (auto& thrd : threads)
        thrd.join();

    return 0;
}