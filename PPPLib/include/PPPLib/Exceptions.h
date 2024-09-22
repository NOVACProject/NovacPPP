#pragma once

#include <exception>
#include <string>

namespace PPPLib
{

class NotFoundException : public std::exception
{
public:
    NotFoundException(std::string msg) : std::exception(), message(msg) {}

    const std::string message;
};


}
