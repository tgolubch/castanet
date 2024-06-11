#!/bin/bash

MyBamPath=$1
echo target_id,target_len && samtools view -H ${MyBamPath} | grep '^@SQ' | cut -d: -f2- | sed s'/\tLN:/,/' | cut -f1,2
