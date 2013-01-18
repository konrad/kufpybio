import xml.etree.ElementTree as ElementTree

class GOXMLParser(object):

    def go_term_information(self, xml):
        tree = ElementTree.fromstring(xml)
        return{
            "name" : tree.findtext("./term/name"),
            "namespace" : tree.findtext("./term/namespace"),
            "def" : tree.findtext("./term/def/defstr")}
