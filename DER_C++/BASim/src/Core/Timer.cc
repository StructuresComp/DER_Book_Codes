#include <iostream>
#include <map>
#include <string>
#include <sys/time.h>
#include <sys/resource.h>
#include "Timer.hh"

namespace BASim {

Timer::TimerMap Timer::timerMap;

inline double timeofday()
{
  timeval tp;
  gettimeofday(&tp, NULL);
  double result = (double) tp.tv_sec + 1e-6 * (double) tp.tv_usec;
  return result;
}

inline double cpuusage()
{
  struct rusage rus;
  getrusage(RUSAGE_SELF, &rus);
  double user = (double) rus.ru_utime.tv_sec + 1e-6
    * (double) rus.ru_utime.tv_usec;
  double sys = (double) rus.ru_stime.tv_sec + 1e-6
    * (double) rus.ru_stime.tv_usec;

  return user + sys;
}

inline double now()
{
  //return timeofday();
  return cpuusage();
}

Timer& Timer::getTimer(std::string name)
{
  TimerMap::iterator iter = timerMap.find(name);
  if (iter != timerMap.end()) return *(iter->second);

  Timer* t = new Timer();
  timerMap[name] = t;
  return *t;
}

void Timer::deleteAllTimers()
{
  TimerMap::iterator iter;

  for (iter = timerMap.begin(); iter != timerMap.end(); ++iter) {
    delete iter->second;
  }

  timerMap.clear();
}

void Timer::report()
{
  TimerMap::iterator iter;

  for (iter = timerMap.begin(); iter != timerMap.end(); ++iter) {
    std::string name = iter->first;
    Timer* t = iter->second;
    std::cout << name << " = " << *t << std::endl;
  }
}

Timer::Timer()
  : startTime(0)
  , elapsed(0)
  , total(0)
  , count(0)
{}

void Timer::clear()
{
  startTime = 0.;
  elapsed = 0.;
}

void Timer::start()
{
  startTime = now();
}

void Timer::stop()
{
  if (startTime) {
    elapsed = now() - startTime;
    total += elapsed;
    times.push_back(elapsed);
    startTime = 0.;
    count++;
  }
}

double Timer::getElapsed()
{
  double e = elapsed;

  if (startTime) {
    e = now() - startTime;
  }

  return e;
}

double Timer::getTotal()
{
  return total;
}

} // namespace BASim
