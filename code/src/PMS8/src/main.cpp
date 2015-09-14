/*
 * main.cpp
 *
 *  Created on: Jun 21, 2012
 *      Author: marius
 */
#include "MotifWorker.h"
#include <ctime>
#include <cmath>
#include "mpi.h"
#include "IndexScheduler.h"

#define MAX_MOTIFS_PER_PROC 100000

struct WorkerStruct {
	IndexScheduler *scheduler;
	MotifWorker *worker;
	int myRank;
	int nTasks;
	clock_t startTime;
};

void doWork(WorkerStruct *ws) {
	int index;
	while ((index = ws->scheduler->requestWork(ws->myRank)) < ws->nTasks) {
		// cout << "Processor " << ws->myRank << " is processing index " << index
		// 		<< endl;
		ws->worker->process(index);
		clock_t endInd = clock();
		long globalTime = (endInd - ws->startTime) / CLOCKS_PER_SEC;
		long projected = (endInd - ws->startTime)
				* ws->nTasks/ (index + 1) / CLOCKS_PER_SEC;
		// cout << "++?! Index " << index << " current time " << (int) globalTime
		// 		<< "s; projected " << projected << "s" << endl;
	}
}

void *runWork(void *ws) {
	doWork((WorkerStruct*) ws);
	return NULL;
}

void createWorkerThread(pthread_t& t, WorkerStruct *ws) {
	pthread_create(&t, NULL, runWork, (void *) ws);
}

int main(int argc, char **argv) {
	int myRank;
	int nProcs; /* number of processes */
	int MotifTransferTag = 2;
	MPI::Status status;

	int totalMsgWords;
	int *b = NULL;

	clock_t startTime = clock();

	MPI::Init_thread(argc, argv, MPI::THREAD_FUNNELED);
	myRank = MPI::COMM_WORLD.Get_rank();
	nProcs = MPI::COMM_WORLD.Get_size();

	cout << "Number of procs " << nProcs << " my id " << myRank << endl;

	if (argc < 4) {
		if (myRank == 0) {
			cout << "Argument: inputFile L d [t [n']]" << endl;
			cout << "inputFile - sequences in fasta format" << endl;
			cout << "L - length of motif" << endl;
			cout
					<< "d - maximum distance between motif and occurence in sequences"
					<< endl;
			cout
					<< "t (optional) - number of lmers in the sets for which we generate common neighborhoods"
					<< endl;
			cout
					<< "n' (optional) - number of sequences for which we initially find motifs"
					<< endl;
		}

		MPI::Finalize();
		exit(1);
	}

	MotifWorker worker;

	if (myRank == 0) {
		int t = -1;
		if (argc > 4)
			t = atoi(argv[4]);
		int nPrime = -1;
		if (argc > 5)
			nPrime = atoi(argv[5]);
		worker.init(atoi(argv[2]), atoi(argv[3]), t, nPrime);
		b = worker.readAndEncodeInput(argv[1], totalMsgWords);
		if (b == NULL) {
			cout << "Unable to open file " << argv[1] << endl;
			MPI::COMM_WORLD.Abort(1);
			exit(1);
		}
	}

	MPI::COMM_WORLD.Bcast(&totalMsgWords, sizeof(int), MPI::CHAR, 0);
	if (myRank != 0) {
//		cout << "Process " << myRank << " has received buffer size "
//				<< totalMsgWords << endl;
		b = new int[totalMsgWords];
	}

	MPI::COMM_WORLD.Bcast(b, totalMsgWords * sizeof(int), MPI::CHAR, 0);
	if (myRank != 0) {
		worker.loadConfigFromBuffer(b);
		cout << "Processor " << myRank << " received ";
		worker.printConfig();
	}

	worker.loadStringsFromBuffer(b);

	int nTasks = worker.getNLmersInFirstString();

	IndexScheduler scheduler(nTasks + nProcs);

	WorkerStruct ws;
	ws.myRank = myRank;
	ws.nTasks = nTasks;
	ws.scheduler = &scheduler;
	ws.worker = &worker;
	ws.startTime = startTime;

	if (myRank == 0 && nProcs > 1) {
		pthread_t thread;
		createWorkerThread(thread, &ws);
		scheduler.loop();
		pthread_join(thread, NULL);
	} else {
		doWork(&ws);
	}

	set<MyString> motifs = worker.getMotifs();
	delete[] b;

	if (myRank == 0) {
		int sz;
		int *buf = worker.allocateMotifBuffer(MAX_MOTIFS_PER_PROC, sz);
		for (int i = 1; i < nProcs; ++i) {
			MPI::COMM_WORLD.Recv(buf, sz, MPI::INT, MPI::ANY_SOURCE,
					MotifTransferTag, status);
			int nm = worker.decodeMotifs(buf, motifs);
			cout << "Processor " << myRank << " received " << nm
					<< " motifs from " << status.Get_source() << endl;
		}
		delete[] buf;

//		vector<MyString> rankedMotifs = worker.rankMotifs(motifs);
		worker.printMotifs(motifs);
		cout << "Total motifs found: " << motifs.size() << endl;
	} else {
		int nMotifs = motifs.size();
		if (nMotifs > MAX_MOTIFS_PER_PROC) {
			cout << "Processor " << myRank << " found " << motifs.size()
					<< " motifs; Keeping first " << MAX_MOTIFS_PER_PROC
					<< " of them" << endl;
			nMotifs = MAX_MOTIFS_PER_PROC;
		}
		int requiredMem;
		int *buf = worker.encodeMotifs(motifs, nMotifs, requiredMem);
		MPI::COMM_WORLD.Send(buf, requiredMem, MPI::INT, 0, MotifTransferTag);
		cout << "Processor " << myRank << " sent " << nMotifs
				<< " motifs to proc 0" << endl;
		delete[] buf;
	}

	clock_t endTime = clock();
	float seconds = (float) (endTime - startTime) / CLOCKS_PER_SEC;
	cout << "Time on processor " << myRank << ": " << seconds << "s; file "
			<< argv[1] << endl;

	/* shut down MPI */
	MPI::Finalize();
	return 0;
}

