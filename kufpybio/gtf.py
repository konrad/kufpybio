from gff3 import Gff3Parser, Gff3Entry

class GtfParser(Gff3Parser):
    
    def _dict_to_entry(self, entry_dict):
        return GtfEntry(entry_dict)

class GtfEntry(Gff3Entry):
    
    def _strip(self, string):
        return string.strip()

    def _attributes(self, attributes_string):
        """Translate the attribute string to dictionary"""

        attributes_string = attributes_string.rstrip(";")
        return dict(
                [self._strip(key_value_pair).split(" ") 
                 for key_value_pair in attributes_string.split(";")])

    def add_attribute(self, key, value):
        self.attributes[key] = value
        self.attribute_string = "; ".join(
            [" ".join(items) for items in self.attributes.items()])
