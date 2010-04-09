// routines for handling terminal access
// provides a terminal class which allows reading of characters
// without stopping the program

#include <stdio.h>
#include <unistd.h>

#include "terminal.hh"

terminal::terminal()
  : _done_init(false)
{
  // we don't change terminal modes if not a terminal
  if( ! isatty( fileno(stdin) ) )
    {
      return;
    }

  FILE* input = fopen("/dev/tty", "r");
  if( input == 0 )
    {
      fprintf(stderr, "[Warning] Could not open /dev/tty\n");
      return;
    }

  // get current attributes
  tcgetattr( fileno(input), &_initial_settings );

  termios new_settings = _initial_settings;

  new_settings.c_lflag &= ~ICANON;
  new_settings.c_lflag &= ~ECHO;
  new_settings.c_cc[VMIN] = 0;
  new_settings.c_cc[VTIME] = 0;
  
  if( tcsetattr( fileno(input), TCSANOW, &new_settings) != 0 )
    {
      fprintf(stderr, "[Warning] Could not set terminal attributes\n");
    }

  fclose(input);

  _done_init = true;
}

terminal::~terminal()
{
  if( ! _done_init )
    return;

  FILE* input = fopen("/dev/tty", "r");
  if( input == 0 )
    {
      fprintf(stderr, "[Warning] Could not open /dev/tty\n");
      return;
    }

  if( tcsetattr( fileno(input), TCSANOW, &_initial_settings) != 0 )
    {
      fprintf(stderr, "[Warning] Could not return terminal attributes\n");
    }

  fclose(input);
}

char terminal::get_char()
{
  if( ! _done_init )
    return 0;

  const int i = fgetc(stdin);
  if( i < 0 )
    return 0;
  else
    return char(i);
}
