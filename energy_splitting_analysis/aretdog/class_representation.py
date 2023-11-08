


class Representation():
    
    # A class that bundles the info on a group and a particular representation of the group, linking characters to conjugacy classes
    # It provides the methods for summing and multiplying representations
    
    def __init__(self, parent_group, rep):
        # initialize either from a dict or an array ordered the same way as the conjugacy classes of parent_group
        self.parent_group = parent_group
        self.characters = {}
        for cc_i in range(len(parent_group.conjugacy_class_names)):
            cc = parent_group.conjugacy_class_names[cc_i]
            if type(rep) == dict:
                self.characters[cc] = rep[cc]
            else:
                # an iterable
                self.characters[cc] = rep[cc_i]
    
    def __str__(self):
        res = f"Rep in group {self.parent_group.name}:"
        for cc in self.characters.keys():
            res += f" '{cc}' = {self.characters[cc]},"
        return(res[:-1])
    def __repr__(self):
        res = f"Rep in group {self.parent_group.name}:"
        for cc in self.characters.keys():
            res += f" '{cc}' = {self.characters[cc]},"
        return(res[:-1])
    
    # ---------------------- rep combination methods
    
    def elementwise_binary_operation(self, other, binary_operation):
        # First, check both reps belong to the same group
        res = {}
        if self.parent_group.name == other.parent_group.name:
            for cc in self.characters.keys():
                res[cc] = binary_operation(self.characters[cc], other.characters[cc])
            return(Representation(self.parent_group, res))
        # Otherwise, we determine if one is in a subgroup of the other. If yes, we return a rep in the subgroup
        if self.parent_group.name in other.parent_group.subgroups:
            other_res_in_self_group = other.parent_group.rep_to_subgroup_rep(self.parent_group.name, other)
            for cc in self.characters.keys():
                res[cc] = binary_operation(self.characters[cc], other_res_in_self_group.characters[cc])
            return(Representation(self.parent_group, res))
        if other.parent_group.name in self.parent_group.subgroups:
            self_res_in_other_group = self.parent_group.rep_to_subgroup_rep(other.parent_group.name, self)
            for cc in other.characters.keys():
                res[cc] = binary_operation(self_res_in_other_group.characters[cc], other.characters[cc])
            return(Representation(other.parent_group, res))
        print("Two provided reps aren't castable to any provided group")
        return(-1)
    
    
    def __add__(self, other):
        return(self.elementwise_binary_operation(other, lambda x, y: x+y))
    def __mul__(self, other):
        return(self.elementwise_binary_operation(other, lambda x, y: x*y))
    def __truediv__(self, other):
        return(self.elementwise_binary_operation(other, lambda x, y: x/y)) #UNSAFE
    def __sub__(self, other):
        return(self.elementwise_binary_operation(other, lambda x, y: x-y))
        

