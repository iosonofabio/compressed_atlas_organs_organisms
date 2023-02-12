# vim: fdm=indent
'''
author:     Fabio Zanini
date:       01/12/22
content:    Load configuration file for a specific instance. This module should
            be imported first to make sure other modules have access to the
            config values.
'''
import yaml

with open("config.yml") as f:
    configuration = yaml.safe_load(f)
