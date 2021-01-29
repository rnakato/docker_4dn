for tag in 18.04 latest #v42 
do
    docker build -t rnakato/4dn:$tag .
    docker push rnakato/4dn:$tag
done
