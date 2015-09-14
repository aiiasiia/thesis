/*
 * IndexScheduler.h
 *
 *  Created on: Nov 24, 2012
 *      Author: marius
 */

#ifndef INDEXSCHEDULER_H_
#define INDEXSCHEDULER_H_
#include "mpi.h"
#include <iostream>
#include <pthread.h>
#include <ctime>
using namespace std;

class IndexScheduler {
public:

	IndexScheduler(int nJobs) :
			WorkRequest(1), nJobs(nJobs), currentIndex(0) {
		pthread_mutex_init(&mutex, NULL);
	}

	~IndexScheduler() {

	}

	int requestWork(int currentThread) {
		if (currentThread == 0) {
			return nextIndex();
		} else {
			MPI::COMM_WORLD.Send(&currentThread, 1, MPI::INT, 0, WorkRequest);
			int index;
			MPI::COMM_WORLD.Recv(&index, 1, MPI::INT, 0, WorkRequest);
			return index;
		}
	}

	void loop() {
		struct timespec tim;
		tim.tv_sec = 0;
		tim.tv_nsec = 10 * 1000000; //10 ms

		MPI::Status status;
		do {

			bool sw = MPI::COMM_WORLD.Iprobe(MPI::ANY_SOURCE, WorkRequest);
			if (sw) {
				int worker;
				MPI::COMM_WORLD.Recv(&worker, 1, MPI::INT, MPI::ANY_SOURCE,
						WorkRequest);

				int msg = nextIndex();
				MPI::COMM_WORLD.Send(&msg, 1, MPI::INT, worker, WorkRequest);
			}

			if (currentIndex >= nJobs)
				return;

			if (!sw) {
//				cout << "Boy am I sleepy" << endl;
				nanosleep(&tim, NULL);
			}
		} while (1);
	}

private:
	const int WorkRequest;
	const int nJobs;
	int currentIndex;
	pthread_mutex_t mutex;

	int nextIndex() {
		int answer;
		pthread_mutex_lock(&mutex);
		answer = currentIndex++;
		pthread_mutex_unlock(&mutex);
		return answer;
	}

}
;

#endif /* INDEXSCHEDULER_H_ */
