# Compressed atlas across organs and organisms
This web application demonstrates the idea of a "compressed cell atlas", i.e. a nonredundant distillation of one or more single cell omics data sets, across organs and organisms.

## Testing
### Virtualenv
To test the web application locally inside a Python virtualenv:

1. make sure you compress the atlases appriopriately or ask Fabio for the compressed files
2. open a terminal on linux/OSX and navigate to the `webapp` subfolder
3. run `./test.sh`

### Docker
To test locally inside a docker container (simulating lightsail virtualization):

1. Test locally as above
2. Start your docker service (e.g. `systemctl start docker`): you probably need superuser rights. If you started it already, you can skip this step.
3. Build (or rebuild) the docker container: `docker build -t compressed-atlas .`
4. Test the image: `docker run -p 5000:5000 compressed-atlas`. In this example it will run the image on port 5000.

## Functionality
Things you can do are shown on the home page once you launch the app, but as a quick (nonexhaustive) summary:
- Show gene expression and chromatin accessibility by cell type
- Show gene expression and chromatin accessibility for a single gene/region across ages in all cell types
- Show gene expression and chromatin accessibility for multiple genes/regions across ages in one cell type
- Look at a disease condition (e.g. lack of exercise)
- Show marker genes/regions for the cell types
- Differential expression/accessibility (WIP)

Unusual features of note:
- The app also tries to give you suggestions for what steps to do next each time you use it.
- Plots can be downloaded in PNG and SVG formats and the data therein as a CSV.
- The compressed atlas itself can be accessed programmatically via the RESTful API, including a full copy of the compressed atlas (WIP)

**NOTE:** raw (uncompressed) data cannot be accessed directly. Ask Emily for access.

## Architecture
The architecture of the compressed atlas is the following:
- A RESTful APIs to request the compressed data (e.g. `/data/gene_names=Col1a1`).
- A set of interative plots (mostly heatmaps or variations on the theme, e.g. dot plots) to visuaise the compressed data.
- A text control system enabling a natural language UX.

At this time, this application is pre-alpha, so the API changes all the time. If you are interested in how it works, write me an email at fabio _DOT_ zanini _AT_ unsw _DOT_ edu _DOT_ au.
