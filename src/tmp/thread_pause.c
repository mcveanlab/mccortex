#include "global.h"
#include "thread_pause.h"
#include <errno.h>

void thread_pause_alloc(ThreadPause *thp)
{
  if(pthread_mutex_init(&thp->pause_lock, NULL) != 0 ||
    pthread_mutex_init(&thp->resume_lock, NULL) != 0 ||
    pthread_mutex_init(&thp->control_lock, NULL) != 0)
  {
    die("pthread_mutex init failed: %s\n", strerror(errno));
  }

  if(pthread_cond_init(&thp->pause_cond, NULL) != 0 ||
     pthread_cond_init(&thp->resume_cond, NULL) != 0)
  {
    die("pthread_cond init failed: %s\n", strerror(errno));
  }

  thp->paused = 0;
  thp->nthreads_running = 0;
}

void thread_pause_dealloc(ThreadPause *thp)
{
  pthread_mutex_destroy(&thp->control_lock);
  pthread_mutex_destroy(&thp->pause_lock);
  pthread_mutex_destroy(&thp->resume_lock);
  pthread_cond_destroy(&thp->pause_cond);
  pthread_cond_destroy(&thp->resume_cond);
}

// Indicate that a thread has started
void thread_pause_started(ThreadPause *thp)
{
  __sync_fetch_and_add(&thp->nthreads_running, 1);
}

// Indicate that a thread has finished
void thread_pause_finished(ThreadPause *thp)
{
  __sync_fetch_and_sub(&thp->nthreads_running, 1);
}

// Returns 1 on success, 0 if someone has already called pause
bool thread_pause_take_control(ThreadPause *thp)
{
  if(pthread_mutex_trylock(&thp->control_lock) != 0) return false;

  thp->paused = true;
  __sync_fetch_and_sub(&thp->nthreads_running, 1);

  pthread_mutex_lock(&thp->pause_lock);
  while(thp->nthreads_running)
    pthread_cond_wait(&thp->pause_cond, &thp->pause_lock);
  pthread_mutex_unlock(&thp->pause_lock);

  return true;
}

// Resume all threads waiting
void thread_pause_release_control(ThreadPause *thp)
{
  thp->paused = false;

  // Wrapping broadcast in lock / unlock required here
  // to avoid:
  // 1: while(thp->paused)
  //                                 2: paused = false; broadcast()
  // 1: wait();
  pthread_mutex_lock(&thp->resume_lock);
  pthread_cond_broadcast(&thp->resume_cond);
  pthread_mutex_unlock(&thp->resume_lock);

  thread_pause_started(thp);
  pthread_mutex_unlock(&thp->control_lock);
}

// Blocks then returns 1 on success, 0 if no one has paused
bool thread_pause_trywait(ThreadPause *thp)
{
  if(!thp->paused) return false;

  __sync_fetch_and_sub(&thp->nthreads_running, 1);

  // Signal wrapped for same reasons as broadcast above
  pthread_mutex_lock(&thp->pause_lock);
  pthread_cond_signal(&thp->pause_cond);
  pthread_mutex_unlock(&thp->pause_lock);

  pthread_mutex_lock(&thp->resume_lock);
  while(thp->paused) pthread_cond_wait(&thp->resume_cond, &thp->resume_lock);
  pthread_mutex_unlock(&thp->resume_lock);

  thread_pause_started(thp);

  return false;
}
