for tag in 20.04 #v42 latest
do
    docker build -t rnakato/4dn:$tag .
    docker push rnakato/4dn:$tag
done
