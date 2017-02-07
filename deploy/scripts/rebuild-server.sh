#!/bin/bash

"$(dirname "$0")"/takedown-serve.sh
"$(dirname "$0")"/images-build.sh
"$(dirname "$0")"/images-push.sh
"$(dirname "$0")"/serve.sh
