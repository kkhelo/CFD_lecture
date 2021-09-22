# include <stdio.h>

int main(void)
{
    int x = 0;
    int y = 0;
    FILE *fp = NULL;
    FILE *gnupipe = NULL;
    char *GnuCommand []= {"set title \"Demo\"", "plot 'data.tmp'"};

    fp = fopen("data.tmp", "w");
    gnupipe = _popen("gnuplot -persitent", "w");

    for(int i = 0; i < 11; i++)
    {
        fprintf(fp, "%d %d\n", x, y);

        x += 1;
        y += 1;
    }

    for (int i = 0; i < 2; i++)
    {
        fprintf(gnupipe, "%s\n", GnuCommand[i]);
    }
}
