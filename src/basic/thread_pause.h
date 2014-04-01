#ifndef THREAD_PAUSE_H_
#define THREAD_PAUSE_H_

// Methods to synchronise a group of threads
// One thread gets control, the reset wait for it to release
// Useful for synchronising threads at checkpoints to e.g. save state

typedef struct
{
  volatile bool paused;
  volatile size_t nthreads_running, nthreads_waiting;
  pthread_mutex_t pause_lock, resume_lock, control_lock;
  pthread_cond_t pause_cond, resume_cond;
} ThreadPause;

void thread_pause_alloc(ThreadPause *thpause);
void thread_pause_dealloc(ThreadPause *thpause);

// Indicate that a thread has started
void thread_pause_started(ThreadPause *thp);
// Indicate that a thread has finished
void thread_pause_finished(ThreadPause *thp);

// Returns 1 on success, 0 if someone has already called
bool thread_pause_take_control(ThreadPause *thpause);
// Resume all threads waiting
void thread_pause_release_control(ThreadPause *thpause);

// Blocks then returns 1, returns 0 immediately if no one has taken control
bool thread_pause_trywait(ThreadPause *thpause);

#endif /* THREAD_PAUSE_H_ */
