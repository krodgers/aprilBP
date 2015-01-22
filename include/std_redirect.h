#include <iostream>
#include <fstream>
#include <sstream>





struct cout_redirect {
    cout_redirect( std::streambuf * new_buffer ) 
        : old( std::cout.rdbuf( new_buffer ) )
    { }

    ~cout_redirect( ) {
        std::cout.rdbuf( old );
    }

private:
    std::streambuf * old;
};

std::stringstream buffer;
std::streambuf *restoreSTDCOUT = std::cout.rdbuf(buffer.rdbuf());

/* std::cout << "Bla" << std::endl; */

/* std::string text = buffer.str(); // text will now contain "Bla\n" */
