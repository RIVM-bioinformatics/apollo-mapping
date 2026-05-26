import argparse

class DynamicHelpTopicAction(argparse.Action):
    def __init__(self, option_strings, dest, **kwargs):
        # Set nargs to 0 because it's a flag that takes no arguments
        super().__init__(option_strings=option_strings, dest=dest, nargs=0, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        """ content is obtained dynamically here """
        content = self.generate_content()
        # if no content, assume generate_content() did the printing itself
        if content:
            print(content)
        # Exit the program after showing the info, behaving like --version or --help
        parser.exit()

    def generate_content(self) -> str:
        """ Define your 'heavy' operation in inheriting classes here

        You can even put heavy imports inside the __call__ method to keep the rest of your script’s global scope clean

        class DynamicHelpTable(DynamicHelpTopicAction):
            generate_content(self):
                return "+---------+---------+\n| Option  | Value   |\n+---------+---------+\n| Dynamic | Data    |\n+---------+---------+"

        parser = argparse.ArgumentParser()
        parser.add_argument(
            '--version',
            action=DynamicHelpTable,
            help="Show dynamic help/version tables and exit"
        )
        """
        raise NotImplementedError

class DynamicHelpTopicShowMarkDownAction(DynamicHelpTopicAction):
    def generate_content(self) -> str:
        # delayed imports ... dunno if preferable
        from rich.console import Console
        from rich.markdown import Markdown
        console = Console()
        md = Markdown(open(self.MARKDOWN_FILE).read())
        console.print(md)
