RUN1='200828_M00777_0134_000000000-D9H23'
RUN2='200829_M00777_0135_000000000-JB82M'

mkdir -p Downloads/miseq_debug/$RUN1
mkdir -p Downloads/miseq_debug/$RUN2

mkdir -p Downloads/miseq_debug/$RUN1/Thumbnail_Images/L001/
mkdir -p Downloads/miseq_debug/$RUN2/Thumbnail_Images/L001/

for entity in InterOp RunInfo.xml runParameters.xml;
do
    echo "copying $entity"
    cp -r NGS/miseq/$RUN1/$entity Downloads/miseq_debug/$RUN1
    cp -r NGS/miseq/$RUN2/$entity Downloads/miseq_debug/$RUN2
done

for cycle in 1 109 110 111 112 113 114 115 116;
do
    echo "copying Thumbnail_Images cycle $cycle"
    cp -r NGS/miseq/$RUN1/Thumbnail_Images/L001/C$cycle.1 Downloads/miseq_debug/$RUN1/Thumbnail_Images/L001/
    cp -r NGS/miseq/$RUN2/Thumbnail_Images/L001/C$cycle.1 Downloads/miseq_debug/$RUN2/Thumbnail_Images/L001/
done

# 1-100
# 101-108
# 109-116
# 116-166