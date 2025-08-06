variable "IMAGE_TAGS" {
  default = "local"
}

function "tag" {
  params = [tag]
  result = "${split(",", "${tag}")}"
}

# Build geo open image
target "lsi" {
  context = "."
  dockerfile = "Dockerfile"
  tags = "${tag("${IMAGE_TAGS}")}"
}