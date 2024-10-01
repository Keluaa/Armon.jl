
barrier:
	spin -a ./barrier.pml
	$(CC) -g -O2 -w -o pan_barrier ./pan.c
	./pan_barrier

ring:
	spin -a ./ring.pml
	$(CC) -g -O2 -w -o pan_ring ./pan.c
	./pan_ring -q

deadlock_too_small:
	spin -a ./deadlock_channel_too_small.pml
	$(CC) -g -O2 -w -o pan_deadlock_channel_too_small ./pan.c
	./pan_deadlock_channel_too_small -q

clean:
	rm -f ./pan* ./*.trail
