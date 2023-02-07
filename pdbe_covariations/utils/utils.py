#!/usr/bin/env python
# software from PDBe: Protein Data Bank in Europe; https://pdbe.org
#
# Copyright 2020 EMBL - European Bioinformatics Institute
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on
# an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied. See the License for the
# specific language governing permissions and limitations
# under the License.

"""Convenience methods
"""

import requests
from bioint.core.exceptions import CovariationsException


def download_file(url, dst):
    """Download file

    Args:
        url (str): What file should be downloaded
        dst (str): Where it should be saved
    """
    response = requests.get(url, allow_redirects=True)

    if response.status_code != 200:
        raise CovariationsException(f"File from {url} could not be downloaded")

    with open(dst, "wb") as fp:
        fp.write(response.content)
