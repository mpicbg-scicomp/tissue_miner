#ifndef TEXT_PROGRESS_BAR
#define TEXT_PROGRESS_BAR

#include <stdio.h>

class TextProgressBar {
public:
  TextProgressBar(const int width=100) : _width(width), _lastPerc(-1), _lastBar(-1) {};

  void update(const double ratio) {
    int perc = int(ratio*100+0.5);
    int bar = int(ratio*_width+0.5);
    if((perc!=_lastPerc) || (bar!=_lastBar)) {
      _lastPerc = perc;
      _lastBar = bar;
      // adapted from http://www.rosshemsley.co.uk:
      printf("%3d%% [", perc);
      for (int x=0; x<bar; x++) printf("#");
      for (int x=bar; x<_width; x++) printf("_");
      // ANSI Control codes to go back to the
      // previous line and clear it.
      // printf("]\n\033[F\033[J");
      // alternative:
      printf("]\r");
      fflush(stdout);
    }
  }
  
  void done(bool remove=true) {
    if(remove) {
      for (int x=0; x<_width+7; x++) printf(" ");
      printf("\r");
      fflush(stdout);
    } else {
      update(1.0);
      printf("\n");
    }
  }
  
private:
  const int _width;
  int _lastPerc, _lastBar;
};

#endif
