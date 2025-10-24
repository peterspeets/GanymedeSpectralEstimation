#include "XML_Interpreter.h"

XML_Interpreter::XML_Interpreter(const string xml_text) {
    /*
    xml_text: a string that is divided into tags.

    The XML_Interpreter object contains the information from a text in the XML format. A hierarchy of tags is
    constructed.
    */
    list<Tag*> activeTags;
    string tag_string;
    int index = 0;
    int i = 0;
    int j = 0;

    while(index < xml_text.length()) {

        if(xml_text[index] == '<') {
            tag_string = Tag::getStringOfFirstTag(xml_text.substr(index,xml_text.length()-index));

            if(tag_string.length() == 0) {
                break;
            }
            Tag* tag = new Tag(tag_string);
            tags.push_back(tag);



            if(! tag->isClosingTag) {
                for (Tag *otherActiveTag : activeTags) {
                    otherActiveTag->innerTags.push_back(tag);
                }
            }
            if(tag->isClosingTag) {
                //for (Tag openingTag : ranges::reverse_view(activeTags)) {
                for (list<Tag*>::reverse_iterator openingTag_iter = activeTags.rbegin(); openingTag_iter  != activeTags.rend(); openingTag_iter ++ ) {
                    Tag* openingTag = *openingTag_iter;
                    tag->label;
                    if(openingTag->label == tag->label && !openingTag->isClosingTag) {
                        openingTag->closingTag = tag;
                        tag->openingTag = openingTag;
                        tag->innerTags = openingTag->innerTags;
                        tag->content = openingTag->content;
                        tag->raw_content = openingTag->raw_content;
                        if(activeTags.back()->label != tag->label) {
                            cout <<"Problem in tags:"<<endl;
                            cout << activeTags.back()->label << endl;
                            cout << tag->label << endl;
                        }

                        activeTags.pop_back();
                        break;

                    }
                }
            }
            for (Tag* otherActiveTag : activeTags) {
                otherActiveTag->raw_content += xml_text.substr(index, tag_string.length());
            }
            if(!tag->isClosingTag) {

                activeTags.push_back(tag);

            }
            index += tag_string.length();


        } else {

            for (Tag* otherActiveTag : activeTags) {

                otherActiveTag->content += xml_text[index];
                otherActiveTag->raw_content += xml_text[index];

            }

            index++;
        }
    }


    for (list<Tag*>::reverse_iterator openingTag_iter = activeTags.rbegin(); openingTag_iter  != activeTags.rend(); openingTag_iter ++ ) {
        Tag* openingTag = *openingTag_iter;
        string xml_tag_string = "</" + openingTag->label + " ";
        for(map<string,string>::iterator it_key = openingTag->keys.begin(); it_key != openingTag->keys.end(); it_key ++) {
            xml_tag_string  += " " + it_key->first + " = \"" + it_key->second  + "\"" ;
        }
        xml_tag_string   += " >";
        Tag* tag = new Tag(xml_tag_string);
        tags.push_back(tag);

    }
    string str_new_content;
    for (Tag* tag : tags) {
        if(tag->label == "Comment") {
            continue;
        }
        str_new_content = "";
        for(i = 0; i < tag->content.length(); i++) {
            if(tag->content[i] == '\\') {
                tag->content[i] = '/';
            }
        }

        for(i = 0; i < tag->content.length(); i++) {
            if(i > 0 && tag->content[i] == '/' && tag->content[i-1] == '/' ) {
                continue;
            } else {
                str_new_content += tag->content[i];
            }
        }
        tag->content = str_new_content ;

    }

}




XML_Interpreter::~XML_Interpreter() {
    for (Tag* tag : tags) {
        delete tag;
    }
    tags.clear();

}

