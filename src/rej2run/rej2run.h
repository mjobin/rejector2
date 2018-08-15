/*
 *  rej2run.h
 *  rej2run
 *
 *  Created by Matt on 12/11/09.
 *  Copyright 2009 Stanford University. All rights reserved.
 *
 */



struct thread_data{
	string outline;
	int xid;
};


class gridthread{
	bool threadalive;
	thread_data tdata;
	pthread_t xthread;
	
public:
	
		//member fxns
	int multithreadinit(thread_data t);
	bool isdead();
	void markdead();
	void cancelthread();
	int whoami();
};



//-------------------
//FUNCTIONS NOT IN CLASSES

void *multi_thread(void *threadarg);