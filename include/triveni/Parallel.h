// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

/// \file Parallel.h
/// \brief Runs independent blocks of work on several threads.
///
/// Each block gets its own random stream (see Random::streamSeed) and writes
/// only its own output slot, so results are identical for any thread count.

#pragma once

#include <algorithm>
#include <atomic>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

namespace triveni
{

/// Calls work(block) for block = 0 .. nBlocks-1 using up to nThreads threads.
/// The first exception thrown by a block is rethrown here.
template <class Work>
void parallelBlocks(int nBlocks, int nThreads, Work&& work)
{
  std::atomic<int> next{0};
  std::exception_ptr failure;
  std::mutex failureMutex;
  auto worker = [&]() {
    for (int block = next++; block < nBlocks; block = next++) {
      try {
        work(block);
      } catch (...) {
        std::lock_guard<std::mutex> lock(failureMutex);
        if (!failure) {
          failure = std::current_exception();
        }
        next = nBlocks; // stop handing out blocks
      }
    }
  };
  const int n = std::clamp(nThreads, 1, std::max(1, nBlocks));
  std::vector<std::thread> threads;
  for (int i = 1; i < n; ++i) {
    threads.emplace_back(worker);
  }
  worker();
  for (auto& t : threads) {
    t.join();
  }
  if (failure) {
    std::rethrow_exception(failure);
  }
}

} // namespace triveni
