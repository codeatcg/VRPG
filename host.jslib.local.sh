
task=$1
function download_file(){
    url=$1
    outfile=$2
    statefile=$3
    if [ ! -f $stateFile ];then
        wget $url -O $outfile
        if [ $? -ne 0 ];then
            echo "Error: file download failed. Please retry later."
            exit 1
        fi
        touch $statefile
    fi
}

if [ "$task" == "local" ];then
    
    if [ -f "templates/vrpg/local.index.html" ];then
        if [ ! -d dependency ];then
            mkdir dependency
        fi
        cd dependency
        
        download_file https://github.com/twbs/bootstrap/releases/download/v4.3.1/bootstrap-4.3.1-dist.zip bootstrap.zip bootstrap.success
        download_file https://github.com/floating-ui/floating-ui/archive/refs/tags/v1.15.0.tar.gz popper.tar.gz popper.success
        download_file https://code.jquery.com/jquery-3.2.1.min.js jquery-3.2.1.min.js jquery.success
        download_file https://raw.githubusercontent.com/lodash/lodash/4.17.10-npm/lodash.min.js lodash.min.js lodash.success
        download_file https://github.com/eligrey/FileSaver.js/archive/refs/tags/v2.0.4.tar.gz FileSaver.tar.gz FileSaver.success
        download_file https://github.com/d3/d3/releases/download/v6.7.0/d3.zip d3.zip d3.success
        download_file https://github.com/tgdwyer/WebCola/archive/refs/tags/v3.3.8.tar.gz WebCola.tar.gz WebCola.success
        
        unzip bootstrap.zip
        tar -zxf popper.tar.gz
        tar -zxf FileSaver.tar.gz
        unzip -d d3 d3.zip
        tar -zxf WebCola.tar.gz
        
        if [ -d ../static/bootstrap-4.3.1-dist ];then
            rm -rf ../static/bootstrap-4.3.1-dist
        fi
        if [ -d ../static/umd ];then
            rm -rf ../static/umd
        fi
        mv bootstrap-4.3.1-dist ../static/bootstrap-4.3.1-dist
        mv floating-ui-1.15.0/dist/umd ../static/umd
        mv jquery-3.2.1.min.js ../static/jquery-3.2.1.min.js
        mv WebCola-3.3.8/WebCola/cola.min.js ../static/cola.min.js
        mv FileSaver.js-2.0.4/dist/FileSaver.js ../static/FileSaver.js
        mv d3/d3.min.js ../static/d3.min.js
        mv lodash.min.js ../static/lodash.min.js

        mv ../templates/vrpg/index.html ../templates/vrpg/cdn.index.html
        mv ../templates/vrpg/local.index.html ../templates/vrpg/index.html
    else
        echo "Nothing was done!"
    fi
elif [ "$task" == "cdn" ];then
    if [ -f "templates/vrpg/cdn.index.html" ];then
        mv templates/vrpg/index.html templates/vrpg/local.index.html
        mv templates/vrpg/cdn.index.html templates/vrpg/index.html
    else
        echo "Nothing was done!"
    fi
else
    echo "Nothing was done!"
fi









