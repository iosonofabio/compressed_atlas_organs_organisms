# vim: fdm=indent
'''
author:     Fabio Zanini
date:       23/04/22
content:    Voice control backend module for compressed atlases
'''
from base64 import decode
import json
from flask import (
        Flask, send_from_directory, request, redirect, url_for,
        Blueprint,
        )
from .interpret_text import interpret_text
from .command_urls import get_command_response


def text_to_response(text_raw):
    '''Convert a text into a JSON response'''

    # Call Google API Speech To Text
    command_dict = interpret_text(text_raw)

    #print('Text raw:', text_raw)
    #print('Command interpreted:', command_dict)

    # Redirect to the correct endpoint
    response = get_command_response(command_dict)

    #print(response)
    return response


mod = Blueprint('text_control_blueprint', __name__)


@mod.route('/submit_text', methods=['POST'])
def text_control():
    text_raw = request.form["text"]
    response = text_to_response(text_raw)

    return response
