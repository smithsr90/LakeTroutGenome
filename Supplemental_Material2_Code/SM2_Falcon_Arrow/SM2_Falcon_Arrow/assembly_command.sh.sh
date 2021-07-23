# Command used to run assembly pipeline
smrtlink-release_6.0.0.47841/smrtcmds/bin/pbsmrtpipe \
pipeline-id pbsmrtpipe.pipelines.polished_falcon_fat \
-e eid_subread:read_metadata.xml \
--preset-json=assembly_settings.json \
--output-dir=assemblyfolder
