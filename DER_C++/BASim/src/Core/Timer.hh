#ifndef TIMER_HH
#define TIMER_HH

#include "BASim/src/Core/Definitions.hh"

namespace BASim {

// macros to compile timing in/out
#ifdef TIMING_ON

#define START_TIMER(name)                         \
  {                                               \
    Timer& t = Timer::getTimer(name); t.start();  \
  }

#define STOP_TIMER(name)                        \
  {                                             \
    Timer& t = Timer::getTimer(name); t.stop(); \
  }                                             \

#define PRINT_TIMER(name)                                         \
  {                                                               \
    std::cout << name << ": " << Timer::getTimer(name).getTotal() \
              << std::endl;                                       \
  }

#define CLEAR_TIMER(name) \
  {                                             \
    Timer::getTimer(name).clear();              \
  }                                             \

#else // timing off
#define START_TIMER(name) {}
#define STOP_TIMER(name) {}
#define PRINT_TIMER(name) {}
#define CLEAR_TIMER(name) {}
#endif

/** A very simple timer class. */
class Timer
{
public:

  typedef std::map<std::string, Timer*> TimerMap;

  /// returns the timer for a given name; reuses timers if they exist
  static Timer& getTimer(std::string name);

  /// deletes all timers created via getTimer() function
  static void deleteAllTimers();

  /// dumps timing information to stdout
  static void report();

  /// reset any previous timings
  void clear();

  /// start the timer
  void start();

  /// stops the timer, adding the elapsed time to the total elapsed time
  void stop();

  /// returns the total elapsed time
  double getElapsed();

  /// returns the total elapsed time
  double getTotal();

  /// returns the number of timings elapsed time represents
  int getCount()
  {
    return count;
  }

  std::vector<double> getTimes()
  {
    return times;
  }

  Timer& operator-(const Timer& t)
  {
    total = total - t.total;
    return *this;
  }

  friend std::ostream& operator<<(std::ostream& os, Timer& timer)
  {
    os << timer.getTotal() << " (" << timer.getCount() << ")";
    return os;
  }

private:
  static TimerMap timerMap; ///< static map of name->Timer

  /// private timer constructor; create timers via static getTimer() function
  Timer();

  std::vector<double> times;

  double startTime;  ///< last time the timer's start function was called
  double elapsed;    ///< time between last start and stop function calls
  double total;      ///< total time between start and stop function calls
  int count;         ///< the number of timings elapsed time represents
};

} // namespace BASim

#endif // TIMER_HH
