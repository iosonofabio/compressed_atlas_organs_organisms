# Create service (only once)
aws lightsail create-container-service --service-name compressed-atlas-service --power small --scale 1

# Check service
aws lightsail get-container-services --service-name compressed-atlas-service

###########################################
# TO UPDATE RUNNING IMAGE, START FROM HERE
###########################################
# Start docker daemon if needed
sudo systemctl start docker
# Make new image (user should be in docker group)
docker build -t compressed-atlas .
# Test the image
docker run -p 5000:5000 compressed-atlas

# Push locally built image to service
aws lightsail push-container-image --service-name compressed-atlas-service --label compressed-atlas --image compressed-atlas
# NOTE: if this is a small change to the app without touching the deps, only the upper
# docker layers should update/be pushed

# Get the new number/name of the image
# e.g. :compressed-atlas-service.compressed-atlas.3

# Change the image referred to by the public containers list
# i.e. edit containers.json with the new image name/number

# NOTE: no need to change the public-endpoint.json unless you want to keep two parallel versions available as backup

# Deploy new image (which is already online and named by now) onto a public URL
aws lightsail create-container-service-deployment --service-name compressed-atlas-service --containers file://containers.json --public-endpoint file://public-endpoint.json

# Check the status of the public container (same command as above)
aws lightsail get-container-services --service-name compressed-atlas-service
# NOTE: the old image version should be RUNNING/ACTIVE, and the service itself should say
# DEPLOYING CONTAINER or somelike like that. After a few minutes, the new version will
# become active

# Check status online (alternative, more visual)
https://lightsail.aws.amazon.com/ls/webapp/us-west-2/container-services/compressed-atlas-service/deployments

# Check out the result (public endpoint)
https://compressed-atlas-service.vphsnedioc8ke.us-west-2.cs.amazonlightsail.com/

# Check the domain
# NOTE: THIS DOES NOT HAVE A DOMAIN YET
#https://neonatallungatlas.org

# Delete service
#aws lightsail delete-container-service --service-name compressed-atlas
