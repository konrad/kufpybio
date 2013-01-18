import os
import urllib.request

class RESTAPI(object):

    """A general class that handles the local file access or the
    retrival of tha file.

    """

    def _get_data(self, path_template, url_template, entity_id):
        file_path = self._file_path(path_template, entity_id)
        if not os.path.exists(file_path):
            self._retrive_data(url_template, entity_id, file_path)
        return(open(file_path).read())

    def _retrive_data(self, url_template, entity_id, file_path):
        data = urllib.request.urlopen(
            self._base_url + url_template % (entity_id)).read()
        data_fh = open(file_path, "wb")
        data_fh.write(data)
        data_fh.close()

    def _file_path(self, path_template, entity_id):
        return(path_template % (self._download_folder, entity_id))

    def _rest_url(self, url_template, entity_id):
        return(url_template % (self._base_url, entity_id))
