#ifndef TAG_H
#define TAG_H
#include <iostream>
#include <string>
#include <list>
#include <map>

using namespace std;

class Tag {
public:
    Tag(const string);
    virtual ~Tag();
    map<string, string> keys;
    list<Tag*> innerTags;
    string label;
    string content;
    string raw_content;
    Tag *closingTag;
    Tag *openingTag;
    static string getStringOfFirstTag(string);


    bool isClosingTag = false;
    bool isOpeningTag = true;

    friend std::ostream& operator<<(std::ostream&, const Tag&);

protected:

private:
};

#endif // TAG_H
