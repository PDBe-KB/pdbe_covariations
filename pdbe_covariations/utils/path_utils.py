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

"""Convenience methods to define file naming scheme across the application
"""


def probability_file(unp_id):
    return f"{unp_id}_prob.txt"


def score_file(unp_id):
    return f"{unp_id}_score.txt"


def probability_log_file(unp_id):
    return f"{unp_id}_prob.log"


def score_log_file(unp_id):
    return f"{unp_id}_score.log"


def msa_file(unp_id):
    return f"{unp_id}.a3m"


def msa_log_file(unp_id):
    return f"{unp_id}_msa.log"


def msa_filtered_file(unp_id, identity, coverage):
    return f"{unp_id}_id{identity}cov{coverage}.a3m"


def msa_filtered_log_file(unp_id):
    return f"{unp_id}_msa_filter.log"


def covariations_file(unp_id):
    return f"{unp_id}_cov.csv"
