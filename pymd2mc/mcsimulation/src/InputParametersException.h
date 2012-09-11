#pragma once

#include <exception>

class InputParametersException : public std::exception
{
    public:
        InputParametersException( std::string reason ) : mReason( reason ) {};
        InputParametersException() : mReason( "" ) {};

        const char * what() const throw()
        {
            return mReason.c_str();
        };
        virtual ~InputParametersException() throw() {};
    protected:
        std::string mReason;
};
