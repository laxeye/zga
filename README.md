# ZGA genome annotator


## Installation

### From Github

Now you can get ZGA by cloning from the repository with `git clone` or by downloading an archive.

### Dependencies

ZGA is written in Python and tested on Python 3.6. ZGA uses several software and libs:

* fastqc
* ea-utils
* bbmap or seqprep
* spades
* unicycler (optional)
* CheckM
* DFast
* BioPython


#### Conda

It's highly recommended to create a new conda environment:

`conda create -n newcoolenv python=3 fastqc ea-utils spades unicycler checkm-genome seqprep dfast bbmap`

and activate it

`conda activate newcoolenv`


Otherwise you may install dependencies to existing conda environment:
`conda install python>=3.6 fastqc ea-utils spades unicycler checkm-genome seqprep dfast bbmap`


#### Other ways

Of course, it's possible to use another ways even compile all tools from source code. In this case you should check if binaries are in your 'PATH' variable.

### Operating systems

ZGA was tested on Ubuntu 18.04. Your feedback is welcome!


## Usage

Run 'zga.py -h' to get a help message.
