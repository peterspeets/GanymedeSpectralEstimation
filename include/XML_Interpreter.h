#ifndef XML_INTERPRETER_H
#define XML_INTERPRETER_H


#include <iostream>
#include <string>
#include <list>
#include "Tag.h"

using namespace std;


class XML_Interpreter {
public:
    XML_Interpreter(const string );
    virtual ~XML_Interpreter();
    list<Tag*> tags;

protected:

private:
};

#endif // XML_INTERPRETER
