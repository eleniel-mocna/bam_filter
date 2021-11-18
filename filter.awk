# Use: awk -v regex=<SEARCHED REGEX> -f filter.awk
# e.g.: awk -v regex=AA[TG]A -f filter.awk 

{   
    if ( $1 ~ "^@") {print $0; next} # this line is a header
    
    if ($1!=old_name)
        {
        if (found==1 && i==2)
            {
            for (k in lanes)
                {
                print lanes[k]
                }
            }
        delete lanes
        i=0
        found=0
        }
    lanes[i]=$0
    i=i+1
    
    if ( $10 ~ regex ) # This line matches given regex        
        {
        found=1        
        }
        
    old_name=$1
}