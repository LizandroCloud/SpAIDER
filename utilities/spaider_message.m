        function [mess] = spaider_message(adv,status)
        
        
                if strcmp(status,'warning')
                    
                ast = repmat('                   *', length(adv),1);
                disp('                   ******************************************************************************');
                disp('                   *                              Warning:                                      *');
                disp('                   *                                                                            *');
                for kj = 1:2;
                    cminus = (76-length(adv{kj'}));
                    message = strcat(ast(kj,:) , adv{kj'}); 
                    disp([message blanks(cminus) '*']); 
                end
                disp('                   *                                                                            *');   
                disp('                   ******************************************************************************'); 
                
                                 
                     mess=0;
                
                
                else
                
                               ast = repmat('                   *', length(adv),1);
                disp('                   ******************************************************************************');
                disp('                   *                              Error:                                        *');
                disp('                   *                                                                            *');
                for kj = 1:2;
                    cminus = (76-length(adv{kj'}));
                    message = strcat(ast(kj,:) , adv{kj'}); 
                    disp([message blanks(cminus) '*']); 
                end
                disp('                   *                                                                            *');   
                disp('                   ******************************************************************************'); 
                
                                 
                     mess=1;
                    
                    
                    
                end
        end