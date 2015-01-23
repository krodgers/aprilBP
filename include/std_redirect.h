#ifdef LOGFILE
#ifndef REDIRECT
#define REDIRECT

#include <iostream>
#include <fstream>
#include <sstream>


/* struct cout_redirect {
    cout_redirect( std::streambuf * new_buffer ) 
        : old( std::cout.rdbuf( new_buffer ) )
    { }

    ~cout_redirect( ) {
        std::cout.rdbuf( old );
    }

private:
    std::streambuf * old;
    };*/

//std::stringstream buffer;
public void redirect_std_out(){
  std::ofstream out("bp_logfile.txt", std::ofstream::app | std::ofstream::out);
  std::streambuf *restoreSTDCOUT = std::cout.rdbuf();
  std::cout.rdbuf(out.rdbuf());
}
/* std::cout << "Bla" << std::endl; */

/* std::string text = buffer.str(); // text will now contain "Bla\n" */
#endif
#endif
