#ifndef TERMINAL__HH
#define TERMINAL__HH

#include <termios.h>

class terminal
{
public:
  terminal();
  ~terminal();

  // return key if pressed or 0 if none
  char get_char();

  // do we have a terminal
  bool is_terminal() const { return _done_init; }

private:
  bool _done_init;
  termios _initial_settings;
};

#endif
