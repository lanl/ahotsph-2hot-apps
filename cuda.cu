/*
 * Copyright 2013 Michael S. Warren.
 * All Rights Reserved.
 */

#include "Msgs.h"
#include "error.h"
#include "mpmy.h"
#include "singlio.h"

extern void initTrace(void);

extern "C" void CUDA_Init(void) {
    int devID;
    struct cudaDeviceProp props;
    cudaError_t err;

    // initTrace();

    devID = 0;
    err = cudaSetDevice(devID);
    if (err != cudaSuccess)
        Error("cudaSetDevice failed, %d %s\n", err, cudaGetErrorString(err));
    err = cudaGetDevice(&devID);
    if (err != cudaSuccess)
        Error("cudaGetDevice failed, %d %s\n", err, cudaGetErrorString(err));
    err = cudaGetDeviceProperties(&props, devID);
    if (err != cudaSuccess)
        Error("cudaGetDeviceProperties failed, %d %s\n", err, cudaGetErrorString(err));
    Msg_do("Compute %d.%d CUDA device %d: [%s]\n", props.major, props.minor, devID, props.name);
    Msg_do("totalGlobalMem:          %.2f GB\n",
           props.totalGlobalMem / (float)(1024 * 1024 * 1024));
    Msg_do("sharedMemPerBlock:       %.2f KB\n", props.sharedMemPerBlock / (float)1024);
    Msg_do("regsPerBlock:            %d\n", props.regsPerBlock);
    Msg_do("warpSize:                %d\n", props.warpSize);
    Msg_do("maxThreadsPerBlock:      %d\n", props.maxThreadsPerBlock);
    Msg_do("totalConstMem:           %.2f KB\n", props.totalConstMem / (float)1024);
    Msg_do("clockRate                %.3f GHz\n", props.clockRate / 1e6);
    Msg_do("memoryClockRate          %.3f GHz\n", props.memoryClockRate / 1e6);
    Msg_do("ECCEnabled               %d\n", props.ECCEnabled);
    Msg_do("computeMode              %d\n", props.computeMode);
    Msg_do("multiProcessorCount      %d\n", props.multiProcessorCount);
    Msg_do("concurrentKernels        %d\n", props.concurrentKernels);
    Msg_do("asyncEngineCount         %d\n", props.asyncEngineCount);
    Msg_do("canMapHostMemory         %d\n", props.canMapHostMemory);
    singlPrintf("Compute %d.%d CUDA device: [%s]\n", props.major, props.minor, props.name);
}
