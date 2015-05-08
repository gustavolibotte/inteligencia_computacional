1 - Baixe o git para Windows pelo endere�o:

https://www.git-scm.com/download/win

2 - Abra o GIT BASH e defina o nome de usuario para os commits:

git config --global user.name "SEU NOME"

3 - Defina tamb�m o e-mail:

git config --global user.email "SEU E-MAIL"

4 - Ainda no GIT BASH, gere uma nova chave SSH:

ssh-keygen -t rsa -C "SEU E-MAIL"

5 - Pressione ENTER para aceitar as op��es:

Enter file in which to save the key (/Users/you/.ssh/id_rsa)
Enter passphrase (empty for no passphrase)
Enter same passphrase again

6 - Ap�s a sua identifica��o ter sido gerada, acione o SSH-AGENT:

ssh-agent -s

7 - Adicione a sua identifica��o ao SSH-AGENT:

ssh-add ~/.ssh/id_rsa

8 - Copie a sua chave p�blica para o clipboard do sistema:

clip < ~/.ssh/id_rsa.pub

9 - Na sua conta do GitHub, acesse as CONFIGURA��ES > SSH keys > Add SSH key,
insira um nome de sua prefer�ncia para a chave e cole-a no campo seguinte.

10 - Abra o GIT BASH e insira:

ssh -T git@github.com

11 - Se a seguinte mensagem for mostrada, a sua chave foi configurada
com sucesso - Responda "yes" para conex�o:

The authenticity of host 'github.com (207.97.227.239)' can't be established.
RSA key fingerprint is 16:27:ac:a5:76:28:2d:36:63:1b:56:4d:eb:df:a6:48.
Are you sure you want to continue connecting (yes/no)? [yes]
Hi username! You've successfully authenticated, but GitHub does not provide shell access.

12 - Com isso, j� estamos prontos para come�ar a controlar o reposit�rio de desenvolvimento.
Para utilizar o git sem complica��es, esse tutorial simplifica muito as coisas:

http://rogerdudler.github.io/git-guide/index.pt_BR.html