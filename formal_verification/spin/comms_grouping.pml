
#define THREADS 4
#define MAX_CYCLE 10


byte start_count = 0;
byte finish_count = 0;
bool can_receive = false;
bool test_lock = false;


proctype remote_process(chan send, recv)
{
    byte data;
    byte cycle = 0;

    do
    :: (cycle == MAX_CYCLE) -> break
    :: (cycle <  MAX_CYCLE) -> {
        recv?data;
        assert(data == cycle);
        data = data + 1;
        cycle = cycle + 1;
        printf("REMOTE: %d to %d\n", data-1, data);
        (can_receive == true) -> send!data;
    }
    od
}


proctype local_thread(chan send, recv; byte tid)
{
    byte data = 0;
    byte recv_data = 0;
    byte cycle = 0;
    byte blocks_done = 0;
    byte blocks_finished = 0;

    do
    :: (cycle == MAX_CYCLE) -> break
    :: (cycle <  MAX_CYCLE) -> {
        /* do some "work" */
        data = cycle;

        /* start_exchange */
        atomic {
            start_count = start_count + 1;
            blocks_done = start_count;
        }
        if
        :: (blocks_done == THREADS) -> {
            send!data;
            printf("LOCAL %d: cycle %d, sent %d\n", tid, cycle, data);
            atomic { start_count = 0; }
        }
        :: else -> skip
        fi
        if
        :: (blocks_done == 1) -> {
            can_receive = true;
        }
        :: else -> skip
        fi

        /* finish_exchange */
        (can_receive && nempty(recv));
        atomic {
            finish_count = finish_count + 1;
            blocks_finished = finish_count;
        }
        if
        :: (blocks_finished == THREADS) -> {
            atomic { finish_count = 0; }
            recv?recv_data;
            printf("LOCAL %d: cycle %d, recv %d\n", tid, cycle, data);
            assert(recv_data == data + 1);
            // can_receive = false;
        }
        :: else -> skip
        fi

        /* end of cycle */
        cycle = cycle + 1
    }
    od
}


init
{
    chan send = [10] of { byte };
    chan recv = [10] of { byte };
    byte i;

    run remote_process(recv, send);
    for (i: 0 .. THREADS-1) {
        run local_thread(send, recv, i);
    }
}

