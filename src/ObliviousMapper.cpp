#include "ObliviousMapper.h"
#include <iostream>
#include <fstream>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/access.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/back_inserter.hpp>

#define MY_THREAD_NUM 48

void client_oblivious_mapper_preprocess(std::vector<uint64_t>& srcPos, std::vector<uint64_t>& dstPos, 
                                    uint64_t plainNumPerPos, uint64_t iter, uint64_t preprocessId, const uint32_t dstTid, bool allowMissing) {
    uint64_t srcPosNum = srcPos.size();
    uint64_t dstPosNum = dstPos.size();
	TaskComm& taskComm = TaskComm::getClientInstance();
	size_t tileNum = taskComm.getTileNum();
	size_t tileIndex = taskComm.getTileIndex();
    CryptoUtil& cryptoUtil = CryptoUtil::getInstance();

    // Build mapping from srcPos value to srcPos index
    std::unordered_map<uint64_t, uint64_t> posMap;
    for (int i=0; i<srcPosNum; ++i) {
        // Insert the first index of a value only
        if (posMap.find(srcPos[i]) == posMap.end())
            posMap.insert({srcPos[i], i});
    }

    taskComm.sendObliviousMapperConfig(srcPosNum, dstPosNum, plainNumPerPos, dstTid);
    
    // Initialize s (random ShareVecVec of length dstPosNum)
    ShareVecVec dstSvv;
    dstSvv.shareVecs.resize(dstPosNum);

    // Recv srcCev
    CipherEntryVec srcCev;
    taskComm.recvCipherEntryVec(srcCev, dstTid);
    // Do mapping & random share splitting
    CipherEntryVec dstCev;
    dstCev.ceVecs.resize(dstPosNum);

    #pragma omp parallel for num_threads(MY_THREAD_NUM)
    for (int i=0; i<dstPosNum; ++i) {
        if (i%10000 == 0) printf("OM preprocessing client split random share (s) i=%d\n", i); 
        auto search = posMap.find(dstPos[i]);
        uint64_t srcIndex = 0;
        if (search == posMap.end()) {
            if (!allowMissing) {
                printf("Did not find pos index in posMap! %lld\n", dstPos[i]);
                exit(-1);
            }
        } else {
            srcIndex = search->second;
        }
        dstCev.ceVecs[i] = srcCev.ceVecs[srcIndex];
        cryptoUtil.splitRandomShareFromCipherEntry(dstCev.ceVecs[i], plainNumPerPos, dstSvv.shareVecs[i].shares);
    }

    // Send dstCev
    taskComm.sendCipherEntryVec(dstCev, dstTid);   

    // Store dstSvv (s)
    std::string fileName = "client_dstSvv_";
    fileName = fileName + std::to_string(iter) + "_" + std::to_string(preprocessId) + "_" + std::to_string(tileIndex) + "_" + std::to_string(dstTid);
    writeShareVecVectoFile(dstSvv, fileName);
}

void server_oblivious_mapper_preprocess(uint64_t iter, uint64_t preprocessId, const uint32_t dstTid) {
    uint64_t srcPosNum = 0;
    uint64_t dstPosNum = 0;
    uint64_t plainNumPerPos = 0;
	TaskComm& taskComm = TaskComm::getServerInstance();
	size_t tileNum = taskComm.getTileNum();
	size_t tileIndex = taskComm.getTileIndex();
    CryptoUtil& cryptoUtil = CryptoUtil::getInstance();

    taskComm.recvObliviousMapperConfig(srcPosNum, dstPosNum, plainNumPerPos, dstTid);

    // Generate r (random ShareVecVec of length srcPosNum)
    ShareVecVec srcSvv;
    srcSvv.shareVecs.resize(srcPosNum);

    #pragma omp parallel for num_threads(MY_THREAD_NUM)
    for (int i=0; i<srcPosNum; ++i) {
        srcSvv.shareVecs[i].shares.resize(plainNumPerPos);
        for (int j=0; j<plainNumPerPos; ++j) {
            srcSvv.shareVecs[i].shares[j] = cryptoUtil.getRandomShare();
        }
    }

    // Encrypt srcSvv
    CipherEntryVec srcCev;
    srcCev.ceVecs.resize(srcPosNum);

    #pragma omp parallel for num_threads(MY_THREAD_NUM)
    for (int i=0; i<srcPosNum; ++i) {
        if (i%10000 == 0) printf("OM preprocessing server encrypt random share ([r]) i=%d\n", i); 
        srcCev.ceVecs[i] = cryptoUtil.encryptToCipherEntry(srcSvv.shareVecs[i].shares, tileIndex);
    }    

    // Send srcCev
    taskComm.sendCipherEntryVec(srcCev, dstTid);

    // Recv dstCev
    CipherEntryVec dstCev;
    taskComm.recvCipherEntryVec(dstCev, dstTid);
    if (dstCev.ceVecs.size() != dstPosNum) {
        printf("Illegal received dstCev size!\n");
        exit(-1);
    }

    // Do decryption
    ShareVecVec dstSvv;
    dstSvv.shareVecs.resize(dstPosNum);

    #pragma omp parallel for num_threads(MY_THREAD_NUM)  
    for (int i=0; i<dstPosNum; ++i) {
        if (i%10000 == 0) printf("OM preprocessing server decrypt random share (M*r-s) i=%d\n", i); 
        cryptoUtil.decryptCipherEntryToInt64Vec(dstCev.ceVecs[i], dstSvv.shareVecs[i].shares);
    }

    // Store dstSvv (M*r-s)
    std::string fileName = "server_dstSvv_";
    fileName = fileName + std::to_string(iter) + "_" + std::to_string(preprocessId) + "_" + std::to_string(tileIndex) + "_" + std::to_string(dstTid);
    writeShareVecVectoFile(dstSvv, fileName);

    // Store srcSvv (r)
    fileName = "server_srcSvv_";
    fileName = fileName + std::to_string(iter) + "_" + std::to_string(preprocessId) + "_" + std::to_string(tileIndex) + "_" + std::to_string(dstTid);
    writeShareVecVectoFile(srcSvv, fileName);
}

uint64_t client_batch_oblivious_mapper_preprocess(std::vector<uint64_t>& srcPos, std::vector<uint64_t>& dstPos, 
                                    uint64_t plainNumPerPos, uint64_t iter, uint64_t preprocessId, const uint32_t dstTid, bool allowMissing) {
    uint64_t srcPosNum = srcPos.size();
    uint64_t dstPosNum = dstPos.size();
	TaskComm& taskComm = TaskComm::getClientInstance();
	size_t tileNum = taskComm.getTileNum();
	size_t tileIndex = taskComm.getTileIndex();
    CryptoUtil& cryptoUtil = CryptoUtil::getInstance();

    uint32_t keySize = cryptoUtil.getKeySize();
    uint32_t plainTextBits = cryptoUtil.getPlainTextBits();
    uint32_t paddingBits = cryptoUtil.getPaddingBits();
    // uint32_t batchSize = keySize / ((plainTextBits + paddingBits) * plainNumPerPos); // Max pos num per cipher entry
    uint32_t batchSize = 5;
    uint32_t plainNumPerCipherEntry = plainNumPerPos * batchSize;

    // Build mapping from srcPos value to srcPos index
    std::unordered_map<uint64_t, uint64_t> posMap;
    for (int i=0; i<srcPosNum; ++i) {
        // Insert the first index of a value only
        if (posMap.find(srcPos[i]) == posMap.end())
            posMap.insert({srcPos[i], i});
    }

    taskComm.sendObliviousMapperConfig(srcPosNum, dstPosNum, plainNumPerPos, dstTid);
    
    // Initialize s (random ShareVecVec of length dstPosNum)
    ShareVecVec dstSvv;
    dstSvv.shareVecs.resize(dstPosNum);

    // Recv srcCev
    CipherEntryVec srcCev;
    taskComm.recvCipherEntryVec(srcCev, dstTid);
    // Do mapping & random share splitting
    CipherEntryVec dstCev;
    dstCev.ceVecs.resize(dstPosNum);

    #pragma omp parallel for num_threads(MY_THREAD_NUM)
    for (int i=0; i<dstPosNum; ++i) {
        if (i%10000 == 0) printf("OM preprocessing client split random share (s) i=%d\n", i); 
        auto search = posMap.find(dstPos[i]);
        uint64_t srcIndex = 0;
        if (search == posMap.end()) {
            if (!allowMissing) {
                printf("Did not find pos index in posMap! %lld\n", dstPos[i]);
                exit(-1);
            }
        } else {
            srcIndex = search->second;
        }
        dstCev.ceVecs[i] = srcCev.ceVecs[srcIndex];
        cryptoUtil.splitRandomShareFromCipherEntry(dstCev.ceVecs[i], plainNumPerCipherEntry, dstSvv.shareVecs[i].shares);
    }

    // Send dstCev
    taskComm.sendCipherEntryVec(dstCev, dstTid);   

    std::cout<<"client dstSvv.shareVecs[0].shares.size() = "<<dstSvv.shareVecs[0].shares.size()<<std::endl;

    for (int i=0; i<batchSize; ++i) {
        // Store dstSvv (s)
        std::string fileName = "client_dstSvv_";
        fileName = fileName + std::to_string(iter + i) + "_" + std::to_string(preprocessId) + "_" + std::to_string(tileIndex) + "_" + std::to_string(dstTid);
        ShareVecVec curDstSvv;
        curDstSvv.shareVecs.resize(dstPosNum);
        for (int j=0; j<dstPosNum; ++j) {
            for (int k=0; k<plainNumPerPos; ++k) {
                curDstSvv.shareVecs[j].shares.push_back(dstSvv.shareVecs[j].shares[i*plainNumPerPos+k]);
            }
        }
        writeShareVecVectoFile(curDstSvv, fileName);
    }

    return batchSize;
}

uint64_t server_batch_oblivious_mapper_preprocess(uint64_t iter, uint64_t preprocessId, const uint32_t dstTid) {
    uint64_t srcPosNum = 0;
    uint64_t dstPosNum = 0;
    uint64_t plainNumPerPos = 0;
	TaskComm& taskComm = TaskComm::getServerInstance();
	size_t tileNum = taskComm.getTileNum();
	size_t tileIndex = taskComm.getTileIndex();
    CryptoUtil& cryptoUtil = CryptoUtil::getInstance();

    taskComm.recvObliviousMapperConfig(srcPosNum, dstPosNum, plainNumPerPos, dstTid);

    uint32_t keySize = cryptoUtil.getKeySize();
    uint32_t plainTextBits = cryptoUtil.getPlainTextBits();
    uint32_t paddingBits = cryptoUtil.getPaddingBits();
    // uint32_t batchSize = keySize / ((plainTextBits + paddingBits) * plainNumPerPos); // Max pos num per cipher entry
    uint32_t batchSize = 5;
    uint32_t plainNumPerCipherEntry = plainNumPerPos * batchSize;

    // Generate r (random ShareVecVec of length srcPosNum)
    ShareVecVec srcSvv;
    srcSvv.shareVecs.resize(srcPosNum);

    #pragma omp parallel for num_threads(MY_THREAD_NUM)
    for (int i=0; i<srcPosNum; ++i) {
        srcSvv.shareVecs[i].shares.resize(plainNumPerCipherEntry);
        for (int j=0; j<plainNumPerCipherEntry; ++j) {
            srcSvv.shareVecs[i].shares[j] = cryptoUtil.getRandomShare();
        }
    }

    // Encrypt srcSvv
    CipherEntryVec srcCev;
    srcCev.ceVecs.resize(srcPosNum);

    #pragma omp parallel for num_threads(MY_THREAD_NUM)
    for (int i=0; i<srcPosNum; ++i) {
        if (i%10000 == 0) printf("OM preprocessing server encrypt random share ([r]) i=%d\n", i); 
        srcCev.ceVecs[i] = cryptoUtil.encryptToCipherEntry(srcSvv.shareVecs[i].shares, tileIndex);
    }    

    // Send srcCev
    taskComm.sendCipherEntryVec(srcCev, dstTid);

    // Recv dstCev
    CipherEntryVec dstCev;
    taskComm.recvCipherEntryVec(dstCev, dstTid);
    if (dstCev.ceVecs.size() != dstPosNum) {
        printf("Illegal received dstCev size!\n");
        exit(-1);
    }

    // Do decryption
    ShareVecVec dstSvv;
    dstSvv.shareVecs.resize(dstPosNum);

    #pragma omp parallel for num_threads(MY_THREAD_NUM)  
    for (int i=0; i<dstPosNum; ++i) {
        if (i%10000 == 0) printf("OM preprocessing server decrypt random share (M*r-s) i=%d\n", i); 
        cryptoUtil.decryptCipherEntryToInt64Vec(dstCev.ceVecs[i], dstSvv.shareVecs[i].shares);
    }

    for (int i=0; i<batchSize; ++i) {
        // Store dstSvv (M*r-s)
        std::string fileName = "server_dstSvv_";
        fileName = fileName + std::to_string(iter + i) + "_" + std::to_string(preprocessId) + "_" + std::to_string(tileIndex) + "_" + std::to_string(dstTid);
        ShareVecVec curDstSvv;
        curDstSvv.shareVecs.resize(dstPosNum);
        for (int j=0; j<dstPosNum; ++j) {
            for (int k=0; k<plainNumPerPos; ++k) {
                curDstSvv.shareVecs[j].shares.push_back(dstSvv.shareVecs[j].shares[i*plainNumPerPos+k]);
            }
        }
        writeShareVecVectoFile(curDstSvv, fileName);

        // Store srcSvv (r)
        fileName = "server_srcSvv_";
        fileName = fileName + std::to_string(iter + i) + "_" + std::to_string(preprocessId) + "_" + std::to_string(tileIndex) + "_" + std::to_string(dstTid);
        ShareVecVec curSrcSvv;
        curSrcSvv.shareVecs.resize(srcPosNum);
        for (int j=0; j<srcPosNum; ++j) {
            for (int k=0; k<plainNumPerPos; ++k) {
                curSrcSvv.shareVecs[j].shares.push_back(srcSvv.shareVecs[j].shares[i*plainNumPerPos+k]);
            }
        }
        writeShareVecVectoFile(curSrcSvv, fileName);
    }

    return batchSize;
}

void client_oblivious_mapper_online(std::vector<uint64_t>& srcPos, std::vector<uint64_t>& dstPos,
                                    const ShareVecVec& srcSvv, ShareVecVec& dstSvv,
                                    uint64_t plainNumPerPos, uint64_t iter, uint64_t preprocessId, const uint32_t dstTid, bool allowMissing) {
    uint64_t srcPosNum = srcPos.size();
    uint64_t dstPosNum = dstPos.size();
	TaskComm& taskComm = TaskComm::getClientInstance();
	size_t tileNum = taskComm.getTileNum();
	size_t tileIndex = taskComm.getTileIndex();
    CryptoUtil& cryptoUtil = CryptoUtil::getInstance();

    // Build mapping from srcPos value to srcPos index
    std::unordered_map<uint64_t, uint64_t> posMap;
    for (int i=0; i<srcPosNum; ++i) {
        // Insert the first index of a value only
        if (posMap.find(srcPos[i]) == posMap.end())
            posMap.insert({srcPos[i], i});
    }

    // Load preprocessed dstSvv (s)
    ShareVecVec preDstSvv;
    std::string fileName = "client_dstSvv_";
    fileName = fileName + std::to_string(iter) + "_" + std::to_string(preprocessId) + "_" + std::to_string(tileIndex) + "_" + std::to_string(dstTid);
    loadShareVecVecFromFile(preDstSvv, fileName);

    // Recv server srcSvv (<x>_1-r)
    ShareVecVec serverSrcSvv;
    taskComm.recvShareVecVec(serverSrcSvv, dstTid);

    // Merge srcSvv with server srcSvv
    ShareVecVec tmpSvv;
    tmpSvv.shareVecs = srcSvv.shareVecs;
    #pragma omp parallel for num_threads(MY_THREAD_NUM)
    for (int i=0; i<srcPosNum; ++i) {
        cryptoUtil.addShareVec(tmpSvv.shareVecs[i].shares, srcSvv.shareVecs[i].shares, serverSrcSvv.shareVecs[i].shares);
    }

    // Do mapping & merge preprocessed dstSvv
    dstSvv.shareVecs.clear();
    dstSvv.shareVecs.resize(dstPosNum);

    #pragma omp parallel for num_threads(MY_THREAD_NUM)
    for (int i=0; i<dstPosNum; ++i) {
        auto search = posMap.find(dstPos[i]);
        uint64_t srcIndex = 0;
        if (search == posMap.end()) {
            if (!allowMissing) {
                printf("Did not find pos index in posMap! %lld\n", dstPos[i]);
                exit(-1);
            }
        } else {
            srcIndex = search->second;
        }
        dstSvv.shareVecs[i].shares.resize(plainNumPerPos);
        cryptoUtil.addShareVec(dstSvv.shareVecs[i].shares, tmpSvv.shareVecs[srcIndex].shares, preDstSvv.shareVecs[i].shares);
    }
}

void server_oblivious_mapper_online(const ShareVecVec& srcSvv, ShareVecVec& dstSvv,
                                    uint64_t iter, uint64_t preprocessId, const uint32_t dstTid) {
    uint64_t srcPosNum = 0;
    uint64_t dstPosNum = 0;
	TaskComm& taskComm = TaskComm::getServerInstance();
	size_t tileNum = taskComm.getTileNum();
	size_t tileIndex = taskComm.getTileIndex();
    CryptoUtil& cryptoUtil = CryptoUtil::getInstance();

    // Load preprocessed srcSvv
    ShareVecVec preSrcSvv;
    std::string fileName = "server_srcSvv_";
    fileName = fileName + std::to_string(iter) + "_" + std::to_string(preprocessId) + "_" + std::to_string(tileIndex) + "_" + std::to_string(dstTid);
    loadShareVecVecFromFile(preSrcSvv, fileName);
    srcPosNum = preSrcSvv.shareVecs.size();

    // Load preprocessed dstSvv ad dstSvv
    fileName = "server_dstSvv_";
    fileName = fileName + std::to_string(iter) + "_" + std::to_string(preprocessId) + "_" + std::to_string(tileIndex) + "_" + std::to_string(dstTid);
    loadShareVecVecFromFile(dstSvv, fileName);
    dstPosNum = dstSvv.shareVecs.size();

    // srcSvv minus preprocessed secSvv
    ShareVecVec tmpSvv;
    tmpSvv.shareVecs = srcSvv.shareVecs;
    #pragma omp parallel for num_threads(MY_THREAD_NUM)
    for (int i=0; i<srcPosNum; ++i) {
        cryptoUtil.substractShareVec(tmpSvv.shareVecs[i].shares, srcSvv.shareVecs[i].shares, preSrcSvv.shareVecs[i].shares);
    }

    // Send server tmpSvv
    taskComm.sendShareVecVec(tmpSvv, dstTid);
}

void writeShareVecVectoFile(ShareVecVec& svv, std::string fileName) {
    TaskComm& taskComm = TaskComm::getClientInstance();
    std::string dirName = "preprocess/";
    dirName += (taskComm.getSetting() + "/");
    fileName = dirName + fileName;
    std::ofstream ofile;
    ofile.open(fileName);
    if (ofile.is_open()) {
        boost::archive::binary_oarchive oa(ofile);  
        oa << svv;
        ofile.flush();
    } else {
        std::cout<<"Unable to open output file "<<fileName<<std::endl;
    }
}

void loadShareVecVecFromFile(ShareVecVec& svv, std::string fileName) {
    TaskComm& taskComm = TaskComm::getClientInstance();
    std::string dirName = "preprocess/";
    dirName += (taskComm.getSetting() + "/");
    fileName = dirName + fileName;
    std::ifstream ifile;
    ifile.open(fileName);
    if (ifile.is_open()) {
        boost::archive::binary_iarchive ia(ifile);
        svv.shareVecs.clear();
        ia >> svv;
    } else {
        std::cout<<"Unable to open input file "<<fileName<<std::endl;
    }
}