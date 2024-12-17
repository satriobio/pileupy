# Installation

## Source

```
git clone https://github.com/satriobio/bulkvis.git 
pip -e .
```

## PIP

Install pileupy using the following command.

```
pip install pileupy
```

## Conda

```
conda install pileupy
```


## Docker

First build the image from the dockerfile.

```
docker build -t pileupy:latest .
```

Run tool from the image

```
docker run -it -v $(pwd):/data/ pileupy:latest
```


