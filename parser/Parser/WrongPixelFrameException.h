/* 
 * File:   WrongPixelFrameException.h
 * Author: mmpi
 *
 * Created on September 25, 2013, 1:39 PM
 */

#ifndef WRONGPIXELFRAMEEXCEPTION_H
#define	WRONGPIXELFRAMEEXCEPTION_H

class WrongPixelFrameException : public std::exception {
  virtual const char* what() const throw() {
    return "Wrong pixel frame!";
  }
};


#endif	/* WRONGPIXELFRAMEEXCEPTION_H */

