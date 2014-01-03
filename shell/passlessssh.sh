ssh-keygen
cat ~/.ssh/id_dsa.pub | ssh bgc@lsr-mat4.nei.nih.gov 'cat >> ~/.ssh/authorized_keys'
