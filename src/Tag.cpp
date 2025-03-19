#include "Tag.h"


using namespace std;



Tag::Tag(const string xmlText) {

    //TODO: enum for these four bools is probably more elegant.
    bool inQuotes = false;
    bool afterEquals = false;
    bool inKeyword = false;
    bool afterKeyword = false;
    string keyword = "";
    string stringInQuotes = "";
    int index = 0;
    int i = 0;
    string xmlString = getStringOfFirstTag(xmlText);

    for(index = 0; index < xmlString.length(); index++) {
        if(xmlString[index] == '>' && !inQuotes) {
            if(label.length() == 0 ) {
                label = keyword;
                keyword = "";
                break;
            }
        }

        if( (xmlString[index] == '<' || xmlString[index] == ' ') && !inQuotes && !inKeyword ) {
            continue;
        }
        if(inKeyword && xmlString[index] == ' ') {
            inKeyword = false;
            afterKeyword = true;
            continue;
        }
        if(xmlString[index] == '\"') {
            if(afterEquals) {
                inQuotes = true;
                afterEquals = false;
            } else if(inQuotes) {
                keys[keyword] = stringInQuotes;
                stringInQuotes = "";
                keyword = "";
                inQuotes = false;

            }
            continue;
        }
        if(inQuotes) {
            stringInQuotes += xmlString[index] ;
            continue;
        }
        if(!inQuotes && xmlString[index] == '=') {
            inKeyword = false;
            afterKeyword = false;
            afterEquals = true;
            continue;
        }
        if(afterKeyword && xmlString[index] != ' ' && xmlString[index] != '=') {
            afterKeyword = false;
            inKeyword = true;
            label = keyword;
            keyword = xmlString[index];
            continue;
        }
        if(!inQuotes && !afterEquals && xmlString[index] != ' ') {
            inKeyword = true;
            keyword += xmlString[index];
            continue;
        }


    }
    if(label[0] ==  '/') {
        isClosingTag = true;
        isOpeningTag = false;
        label = label.substr(1,label.length()-1);
    } else {
        isClosingTag = false;
        isOpeningTag = true;

    }

    bool selfClosing = false;
    for(i = label.size()-1; i > 0 ; i--) {
        if(xmlString[i] == ' ' || xmlString[i] == '>'  ) {
            continue;
        } else if(xmlString[i] == '/') {
            selfClosing = true;
            break;
        } else {
            break;
        }

    }

    if(label[label.size()-1] ==  '/') {
        selfClosing = true;
    }

    if( selfClosing) {
        isClosingTag = true;
        isOpeningTag = true;
    }


}

string Tag::getStringOfFirstTag(const string xml_string) {
    bool inQuotes = false;
    int i;
    int begin_index = 0;
    int end_index = 0;
    for(i = 0; i < xml_string.size(); i++ ) {
        if(xml_string[i] == '<') {
            begin_index = i;
            break;
        }
    }

    for(i = 0; i < xml_string.size(); i++ ) {
        if(xml_string[i] ==  '\"') {
            inQuotes = !inQuotes;
        }

        if(xml_string[i] == '>') {
            end_index = i;
            break;
        }
    }
    return xml_string.substr(begin_index,end_index - begin_index+1);
}

ostream& operator<<(ostream& os, const Tag& self) {
    os << "<";
    if(self.isClosingTag) {
        os << "/";
    }
    os << self.label << " " ;
    for (const pair<string, string>& pair : self.keys) {
        os << pair.first << " = \"" << pair.second << "\" ";
    }
    os << ">";

    return os;
}

Tag::~Tag() {
    //dtor
}

