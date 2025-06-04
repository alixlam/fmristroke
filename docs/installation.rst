
.. _Installation:

-------------
Installation
-------------


Manually Prepared Environment (Python 3.9+)
===========================================

On a functional Python 3.9 (or above but less than 3.12) environment with ``pip`` installed,
*fMRIStroke* can be installed using the habitual command ::

    $ python -m pip install git+https://github.com/alixlam/fmristroke.git 



External Dependencies
---------------------
*fMRIStroke* is written using Python 3.9 (or above but less than 3.12), and is based on
nipype_.

Containerized execution (Docker and Singularity)
================================================
Docker and Singularity are containerization technologies that allow you to run applications in isolated environments, ensuring that all dependencies are met without affecting the host system.

Docker
------
For every new version that is released, a corresponding Docker image is generated. The Docker image becomes a container when the execution engine loads the image and adds an extra layer that makes it runnable. In order to run the Docker images, the Docker Engine must be installed.

1. **Install Docker**

   If not already installed, refer to the `official Docker installation guide <https://docs.docker.com/get-docker/>`_.

2. **Pull the Docker image** ::

       $ docker pull alixlam/fmristroke:<latest-version>

    Replace ``<latest-version>`` with the desired tag (e.g., ``latest``).

3. **Run the Docker container** ::
         $ docker run -ti --rm \
                -v path/to/data:/data:ro \        # read-only, for data
                -v path/to/output:/out \          # read-write, for outputs
                alixlam/fmristroke:<latest-version> \
                /data /out/out \
                participant


Singularity
---------
Singularity is an alternative to Docker, particularly suited for high-performance computing (HPC) environments where users may not have root access. To run the *fMRIStroke* container with Singularity, follow these steps:

1. **Install Singularity**

   If not already installed, refer to the `official Singularity installation guide <https://docs.sylabs.io/guides/3.5/user-guide/introduction.html>`_.

2. **Build the Singularity image from Docker** ::

       $ singularity build fmristroke.sif docker://alixlam/fmristroke:<latest-version>

   Replace ``<latest-version>`` with the desired tag (e.g., ``latest``).

3. **Run the Singularity container** ::

       $ singularity run \
           -B /path/to/data:/data:ro \
           -B /path/to/output:/out \
           fmristroke.sif \
           /data /out/out \
           participant

   - ``-B`` binds directories from the host system into the container.
   - ``/data`` and ``/out`` refer to the internal container paths.

4. **Notes**:
   - Singularity does not require root privileges to execute containers.
   - It integrates well with shared filesystems and job schedulers.
   - You may use ``singularity exec`` instead of ``run`` for more control over the execution command.
